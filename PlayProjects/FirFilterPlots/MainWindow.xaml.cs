using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Media;

using MathNet.Filtering.FIR;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.Providers.FourierTransform;

using Plot2D_Embedded;

namespace FirFilterPlots
{
    public partial class MainWindow : System.Windows.Window
    {
        //*******************************************************************

        // specify input time domain signal
        private readonly double       inputSampleRate = 4096; 
        private readonly List<double> frequencies = new List<double> () {500};
        private readonly double       plotDuration = 0.05;

        //*******************************************************************

        // store input time domain signal
        private List<Point> inputSignal = new List<Point> (); // (sampleTime, sampleAmplitude)

        //*******************************************************************

        // specify filter & FFT
        private readonly double cutoffFrequency = 1024; // Hz
        private readonly int decimation = 2;
        private readonly int fftSize = 1024;

        //*******************************************************************

        // store the filter
        private double [] filterCoefs; 
        private OnlineFirFilter filter;

        private List<Point> filterSamples;
        private readonly List<Point> filterSpectrum = new List<Point> ();

        //*******************************************************************

        // store the filtered time domain signal
        private List<Point> filteredSignal; // = new List<Point> (); // (sampleTime, sampleAmplitude)

        //*******************************************************************

        // FFT results
        private readonly List<Point> inputSpectrum    = new List<Point> (); // (frequency, magnitude)
        private readonly List<Point> filteredSpectrum = new List<Point> ();

        //*******************************************************************
        //*******************************************************************
        //*******************************************************************

        public MainWindow ()
        {
            InitializeComponent ();
            timeDomainFigure.PlotWindowReady += PlotReady;
            freqDomainFigure.PlotWindowReady += PlotReady;
        }

        int readyCntr = 0;

        private void PlotReady (object sender)
        {
            readyCntr++;

            if (readyCntr == 2)
            { 
                DoCalculations ();
                TimePlotSelectCB_Click (null, null); // draw the default plot selections
                FreqPlotSelectCB_Click (null, null);
            }

            if (sender is Bare2DPlot fig)
                fig.AxesEqual = false;
        }

        //****************************************************************************************

        private void DoCalculations ()
        {
            try
            {
                filterCoefs   = CalculateFilterCoefficients (inputSampleRate, cutoffFrequency);
              //filterCoefs   = GetFilterCoefficients ();

                filter        = GenerateFilter (filterCoefs);
                filterSamples = PadFilterCoeffs (inputSampleRate, filterCoefs, fftSize);


                inputSignal    = GenerateInputSignal (fftSize, filterCoefs.Length, decimation);
                filteredSignal = RunFilter (inputSignal, filter, decimation);
                CalculateSpectra ();

                PlotInputSignal ();

                Console.WriteLine (filterCoefs.Length + " filter coefficients");
                Console.WriteLine ("filtered signal count = " + filteredSignal.Count);
        
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception in DoCalculations: " + ex.Message);
              //Console.WriteLine ("Exception in DoCalculations: " + ex.StackTrace);
            }
        }

        //***************************************************************************************
        //
        // make input signal, filter coefficients and output signal
        //
        private List<Point> GenerateInputSignal (int fftSize, int filterLength, int filterDecimation)
        {
            List<Point> inputSignal = new List<Point> ();
            Random rand = new Random ();

            int sampleCount = fftSize * filterDecimation + filterLength;
            Console.WriteLine ("SampleCount = " + sampleCount.ToString ());

            double t = 0;

            for (int i = 0; i<sampleCount; i++, t+=1/inputSampleRate)
            {
                double s = 0;
                foreach (double f in frequencies)
                    s += Math.Sin (2 * Math.PI * f * t) + (rand.NextDouble () - 0.5) * 0.0001;
                
                inputSignal.Add (new Point (t, s));
            }

            return inputSignal;
        }

        //***************************************************************************************
        //
        // GenerateFilterCoefs ()
        //

        OnlineFirFilter GenerateFilter (double [] filterCoefs)
        {
            return new OnlineFirFilter (filterCoefs);
        }

        List<Point> PadFilterCoeffs (double sampleRate, double [] filterCoefs, int fftSize)
        {
            List<Point> filterSamples = new List<Point> ();
            int i;

            for (i=0; i<filterCoefs.Length; i++)
                filterSamples.Add (new Point (i / sampleRate, filterCoefs [i]));

            for ( ; i<fftSize; i++)
                filterSamples.Add (new Point (i / sampleRate, 0));

            return filterSamples;
        }

        #pragma warning disable IDE0051 // Remove unused private members

        double [] CalculateFilterCoefficients (double sampleRate, double cutoffFreq)
        {
            return FirCoefficients.LowPass (sampleRate, cutoffFrequency); 
        }

        double [] GetFilterCoefficients () // http://t-filter.engineerjs.com/
        {
            return new double [] 
                                 //{-0.02010411882885732,
                                  // -0.05842798004352509,
                                  // -0.061178403647821976,
                                  // -0.010939393385338943,
                                  //  0.05125096443534972,
                                  //  0.033220867678947885,
                                  // -0.05655276971833928,
                                  // -0.08565500737264514,
                                  //  0.0633795996605449,
                                  //  0.31085440365663597,
                                  //  0.4344309124179415,
                                  //  0.31085440365663597,
                                  //  0.0633795996605449,
                                  // -0.0856550073726451,
                                  // -0.05655276971833928,
                                  //  0.033220867678947885,
                                  //  0.05125096443534972,
                                  // -0.010939393385338943,
                                  // -0.061178403647821976,
                                  // -0.05842798004352509,
                                  // -0.02010411882885732};

                                  { 0.020132210722515607,
                                    0.014337588088026872,
                                    -0.06042518986827016,
                                    -0.11688176581198412,
                                    0.015390548525687591,
                                    0.30600043556088063,
                                    0.464289723357815,
                                    0.30600043556088063,
                                    0.015390548525687591,
                                    -0.11688176581198412,
                                    -0.06042518986827016,
                                    0.014337588088026872,
                                    0.020132210722515607 };


        }

        #pragma warning restore IDE0051 // Remove unused private members

        //***************************************************************************************
        //
        // RunFilter
        //

        /// <summary>
        /// Wrapper for MathNet FIR filter ProcessSamples
        /// </summary>
        /// <param name="inputSignal"></param>
        /// <param name="outputSignal"></param>
        /// <param name="filter"></param>
        /// <param name="decimation"></param>
        /// 

        List<Point> RunFilter (List<Point> inputSignal, OnlineFirFilter filter, int decimation)
        {
            List<Point> outputSignal = new List<Point> ();

            // extract samples without time tag
            double [] samples = new double [inputSignal.Count];

            for (int i=0; i<inputSignal.Count; i++)
                samples [i] = inputSignal [i].Y;

            // fullOutputSignal contains startup transients and has not been decimated
            double [] fullOutputSignal = filter.ProcessSamples (samples);

            // remove startup transients
            double [] trimmedOutputSignal = new double [fullOutputSignal.Length - filterCoefs.Length];

            for (int i=0; i<trimmedOutputSignal.Length; i++)
                trimmedOutputSignal [i] = fullOutputSignal [i + filterCoefs.Length];

            // decimate and convert to plot format
            for (int i = 0; i<trimmedOutputSignal.Length; i+=decimation)
                outputSignal.Add (new Point (inputSignal [i].X, trimmedOutputSignal [i]));

            return outputSignal;
        }

        //******************************************************************

        void CalculateSpectra ()
        {
            Console.WriteLine ("input FFT Resolution = " + inputSampleRate / fftSize);
            Console.WriteLine ("filtered FFT Resolution = " + inputSampleRate / decimation / fftSize);

            RunFFT (inputSignal,    inputSpectrum);
            RunFFT (filteredSignal, filteredSpectrum);
            RunFFT (filterSamples,  filterSpectrum);
        }


        private void PlotInputSpectrum ()
        { 
            LineView lv = new LineView (inputSpectrum);
            lv.Color = Brushes.Red;
            freqDomainFigure.Plot (lv);
            freqDomainFigure.RectangularGridOn = true;
        }

        private void PlotOutputSpectrum ()
        { 
            LineView lv2 = new LineView (filteredSpectrum);
            lv2.Color = Brushes.Green;
            freqDomainFigure.Plot (lv2);
            freqDomainFigure.RectangularGridOn = true;
        }

        private void PlotFilterSpectrum ()
        { 
            LineView lv3 = new LineView (filterSpectrum);
            lv3.Color = Brushes.Blue;
            freqDomainFigure.Plot (lv3);
            freqDomainFigure.RectangularGridOn = true;
        }

        //**************************************************************************
        //
        // RunFFT ()
        //

        /// <summary>
        ///  Wrapper for MathNet FFT. Input Points are (time, ampl), Output Points are (freq, magnitude)
        /// </summary>
        /// <param name="timeSignal"></param>
        /// <param name="magSpectum"></param>
        /// <param name="HammWindow"></param>
        
        void RunFFT (List<Point> timeSignal, List<Point> magSpectum, bool HammWindow = false)
        {
            try
            { 
                int sampleCount = fftSize;
                int pad = sampleCount.IsEven () ? 2 : 1;
                int workBufferSize = sampleCount + pad;

                double [] workBuffer = new double [workBufferSize]; // before FFT: input samples
                                                                    // after FFT: half of complex spectrum
                int get = timeSignal.Count - fftSize;

                // copy input to WorkBuffer. the "pad" words were initialized to 0 and will remain 0
                int ii;

                for (ii=0; ii<sampleCount; ii++)
                    workBuffer [ii] = timeSignal [get + ii].Y;

                // optional windowing
                if (HammWindow == true)
                {
                    double [] Hamm = MathNet.Numerics.Window.Hamming (sampleCount);

                    for (int i = 0; i<fftSize; i++)
                        workBuffer [i] *= Hamm [i];
                }

                Fourier.ForwardReal (workBuffer, fftSize, FourierOptions.Matlab); // .NoScaling);

                double [] fftReal = new double [fftSize]; // complex spectrum unpacked and reflected into these
                double [] fftImag = new double [fftSize];

                int put = 0;

                for (int k = 0; k<workBuffer.Length; k+=2, put++)
                {
                    fftReal [put] = workBuffer [k];
                    fftImag [put] = workBuffer [k+1];
                }

                put = fftReal.Length - 1;

                for (int k = 2; k<workBuffer.Length; k+=2, put--)
                {
                    fftReal [put] = workBuffer [k];
                    fftImag [put] = workBuffer [k+1] * -1;
                }

                double sampleFreq = 1 / (timeSignal [1].X - timeSignal [0].X);

                double [] frequencyScale = Fourier.FrequencyScale (fftSize, sampleFreq);
                double [] magnitude = new double [fftSize];

                for (int i = 0; i<fftSize; i++)
                    magnitude [i] = 20 * Math.Log10 (Math.Sqrt (fftReal [i] * fftReal [i] + fftImag [i] * fftImag [i]));

                int L2 = 1 + fftSize / 2;

                for (int i = L2; i<fftSize; i++)
                    magSpectum.Add (new Point (frequencyScale [i], magnitude [i]));

                for (int i = 0; i<L2; i++)
                    magSpectum.Add (new Point (frequencyScale [i], magnitude [i]));
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception in RunFFT: " + ex.Message);
            }
        }

        //************************************************************************
        //
        // Plot signals
        //
        private void PlotInputSignal ()
        {
            LineView lv1 = new LineView (inputSignal);
            lv1.Color = Brushes.Red;
            timeDomainFigure.Plot (lv1);

            // limit display for better view of first part
            timeDomainFigure.GetAxes (out _, out _, out double ymin, out double ymax);
            timeDomainFigure.SetAxes (0, plotDuration, ymin, ymax);

            timeDomainFigure.RectangularGridOn = true;
        }

        private void PlotOutputSignal ()
        {
            LineView lv2 = new LineView (filteredSignal);
            lv2.Color = Brushes.Green;
            timeDomainFigure.Plot (lv2);

            // limit display for better view of first part
            timeDomainFigure.GetAxes (out _, out _, out double ymin, out double ymax);
            timeDomainFigure.SetAxes (0, plotDuration, ymin, ymax);

            timeDomainFigure.RectangularGridOn = true;
        }

        //************************************************************************
        //
        // Plot filter coefficients
        //
        private void PlotFilterCoefficients ()
        {
            List<Point> pts = new List<Point> ();

            for (int i = 0; i<filterCoefs.Length; i++)
                pts.Add (new Point (i, filterCoefs [i]));

            LineView lv = new LineView (pts);

            timeDomainFigure.Plot (lv);
            timeDomainFigure.RectangularGridOn = true;
        }

        //*****************************************************************************

        private void TimePlotSelectCB_Click (object sender, RoutedEventArgs e)
        {
            timeDomainFigure.Clear ();

            if (PlotInputCB.IsChecked  == true) PlotInputSignal ();
            if (PlotOutputCB.IsChecked == true) PlotOutputSignal ();
            if (PlotFilterCB.IsChecked == true) PlotFilterCoefficients ();

            timeDomainFigure.RectangularGridOn = true;
        }

        private void FreqPlotSelectCB_Click (object sender, RoutedEventArgs e)
        {
            freqDomainFigure.Clear ();

            if (FreqPlotInputCB.IsChecked  == true) PlotInputSpectrum ();
            if (FreqPlotOutputCB.IsChecked == true) PlotOutputSpectrum ();
            if (FreqPlotFilterCB.IsChecked == true) PlotFilterSpectrum ();

            freqDomainFigure.RectangularGridOn = true;
        }

        //*****************************************************************************

        private void ZoomOptionButton_Checked (Bare2DPlot fig, RadioButton rb)
        {
            string tag = rb.Tag as string;

            if (tag == null)
                throw new Exception ("Zoom radio button tag read failed");

            switch (tag)
            {
                case "Zoom_Both": fig.ZoomBoth = true; break;
                case "Zoom_X":    fig.ZoomXOnly = true; break;
                case "Zoom_Y":    fig.ZoomYOnly = true; break;
                default:          throw new Exception ("Invalid zoom option");
            }
        }

        private void TimeZoomOptionButton_Checked (object sender, RoutedEventArgs args)
        {
            ZoomOptionButton_Checked (timeDomainFigure, sender as RadioButton);
        }

        private void FreqZoomOptionButton_Checked (object sender, RoutedEventArgs args)
        {
            ZoomOptionButton_Checked (freqDomainFigure, sender as RadioButton);
        }

        //**************************************************************************

        private void Window_Loaded (object sender, RoutedEventArgs e)
        {
            TimeZoomX_Button.IsChecked = true;
            FreqZoomX_Button.IsChecked = true;
        }
    }
}