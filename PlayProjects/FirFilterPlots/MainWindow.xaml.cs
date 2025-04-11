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
    #pragma warning disable IDE0051 // Remove unused private members

    public partial class MainWindow : System.Windows.Window
    {
        //*******************************************************************

        // specify input time domain signal
        private readonly double       inputSampleRate = 100000; 
        //private readonly List<double> frequencies = new List<double> () {1000, 5000, 10000, 15000, 20000, 25000 };
        private readonly List<double> frequencies = new List<double> () {100, 19600 };
        private readonly double       plotDuration = 0.05;

        //*******************************************************************

        // store input time domain signal
        private List<Point> inputSignal = new List<Point> (); // (sampleTime, sampleAmplitude)

        //*******************************************************************

        // specify filter & FFT
        private readonly double cutoffFrequency = 13000; //3125;// 50000 / 16; // Hz
        private readonly int decimation = 1; // 16;
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

                //***********************************************************

                double sum = 0;

                foreach (double d in filterCoefs)
                    sum += d;
        
                Console.WriteLine ("filter coef sum = " + sum.ToString ());

                for (int i=0; i<filterCoefs.Length; i++)
                    filterCoefs [i] /= sum;

                //***********************************************************

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
                    s += Math.Cos (2 * Math.PI * f * t) + (rand.NextDouble () - 0.5) * 0.0001;
                
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

        double [] CalculateFilterCoefficients (double sampleRate, double cutoffFreq)
        {
            return FirCoefficients.LowPass (sampleRate, cutoffFrequency); 
        }

        double [] GetFilterCoefficients () // http://t-filter.engineerjs.com/
        {
            return new double []
                                 {
                                    -0.007193104140960058,
                                    -0.005690941419964349,
                                    -0.006732491431624812,
                                    -0.006734559847354272,
                                    -0.005157607255598935,
                                    -0.0015142487868529896,
                                    0.0045471213869874285,
                                    0.013181429570356848,
                                    0.02427973460356201,
                                    0.037423914508760085,
                                    0.05195106974133303,
                                    0.0669359545576894,
                                    0.08131554113691528,
                                    0.09397473286411374,
                                    0.10387823859496062,
                                    0.11018591835898431,
                                    0.11234988206563754,
                                    0.11018591835898431,
                                    0.10387823859496062,
                                    0.09397473286411374,
                                    0.08131554113691528,
                                    0.0669359545576894,
                                    0.05195106974133303,
                                    0.037423914508760085,
                                    0.02427973460356201,
                                    0.013181429570356848,
                                    0.0045471213869874285,
                                    -0.0015142487868529896,
                                    -0.005157607255598935,
                                    -0.006734559847354272,
                                    -0.006732491431624812,
                                    -0.005690941419964349,
                                    -0.007193104140960058
                                   };
        }

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
            Console.WriteLine ("input FFT Resolution = " + inputSampleRate / fftSize + " Hz / bin");
            Console.WriteLine ("filtered FFT Resolution = " + inputSampleRate / decimation / fftSize + " Hz / bin");

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

        private void PlotTimeSignal (List<Point> signal, Brush color)
        { 
            //List<double> amplitudeOnly = new List<double> ();

            //foreach (Point pt in signal)
            //{ 
            //    amplitudeOnly.Add (pt.Y);

            //    if (amplitudeOnly.Count > 200)
            //        break;
            //}




          LineView lv1 = new LineView (signal);
          //  LineView lv1 = new LineView (amplitudeOnly);

            lv1.Color = color;
            timeDomainFigure.Plot (lv1);

            // limit display for better view of first part
            timeDomainFigure.GetAxes (out _, out _, out double ymin, out double ymax);
            timeDomainFigure.SetAxes (0, plotDuration, ymin, ymax);

            timeDomainFigure.RectangularGridOn = true;
        }

        private void PlotInputSignal ()
        {
            PlotTimeSignal (inputSignal, Brushes.Red);
        }

        private void PlotOutputSignal ()
        {
            PlotTimeSignal (filteredSignal, Brushes.Green);
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

    #pragma warning restore IDE0051 // Remove unused private members
}