using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Media;

using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using Plot2D_Embedded;

namespace RealFFT
{
    public partial class MainWindow : System.Windows.Window
    {
        const int sampleCount = 32; // 1024;

        const double sampleFreq = 2048; // samples per second
        const double resolution = sampleFreq / sampleCount;

        double [] inputSignal = new double [sampleCount]; // real time samples
        double [] sampleTime  = new double [sampleCount];

        double [] fftReal     = new double [sampleCount]; // complex spectrum unpacked and reflected into these
        double [] fftImag     = new double [sampleCount];

        List<double> frequencies = new List<double> () {0, 4 * resolution};
        List<double> amplitudes  = new List<double> () {1, 1, 0.5}; // one for each frequency

        public MainWindow ()
        {
            InitializeComponent ();
            figure.PlotWindowReady += PlotFigureReady;
        }

        private void PlotFigureReady (object sender)
        {
            DoCalculations ();
            PlotSelectCB_Click (null, null); // draw the default plot selections
        }

        //*****************************************************************************

        private void DoCalculations ()
        {
            try
            { 
                //
                // sample times
                //
                for (int i = 0; i<sampleCount; i++)
                    sampleTime [i] = i / sampleFreq;

                //
                // real input signal
                //
                Console.WriteLine ("sample freq = {0:0.###} sa/sec", sampleFreq);
                Console.WriteLine ("resolution = {0:0.###} Hz", resolution);

                Console.Write ("Signal present in bins: ");
                for (int j = 0; j<frequencies.Count; j++)
                    Console.Write ("{0:0.#}, ", frequencies [j] / resolution);
                Console.WriteLine ('\n');

                for (int i = 0; i<sampleCount; i++)
                {
                    inputSignal [i] = 0;

                    for (int j = 0; j<frequencies.Count; j++)
                        inputSignal [i] += amplitudes [j] * Math.Cos (2 * Math.PI * frequencies [j] * sampleTime [i] + 10 * Math.PI / 180);
                }

                //
                // run FFT
                //
                int pad = sampleCount.IsEven () ? 2 : 1;

                double [] workBuffer = new double [sampleCount + pad]; // before FFT: input signal
                                                                       // after FFT: half of complex spectrum
                inputSignal.CopyTo (workBuffer, 0);

                Fourier.ForwardReal (workBuffer, sampleCount, FourierOptions.NoScaling);

                int put = 0;
                
                for (int k=0; k<workBuffer.Length; k+=2, put++)
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



                for (int k=0; k<sampleCount; k++)
                    Console.WriteLine ("{0:0}: {1:0.000}, {2:0.000}", k, fftReal [k], fftImag [k]);
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception in DoCalculations: " + ex.Message);
            }
        }

        //*****************************************************************************

        private void PlotSelectCB_Click (object sender, RoutedEventArgs e)
        {
            figure.Clear ();

            if (PlotSignalCB.IsChecked == true) PlotSignal ();
            if (PlotSpectrumCB.IsChecked == true) PlotSpectrum ();
            //if (PlotFilterCB.IsChecked == true) PlotFilterCoefficients ();

          //  figure.RectangularGridOn = true;
        }

        //*****************************************************************************

        // swap halves

        private void PlotSpectrum ()
        {
            try
            { 
                List<Point> pts = new List<Point> ();

                double [] frequencyScale = Fourier.FrequencyScale (sampleCount, sampleFreq);
                double [] magnitude = new double [sampleCount];

                for (int i = 0; i<sampleCount; i++)
                    magnitude [i] = Math.Sqrt (fftReal [i] * fftReal [i] + fftImag [i] * fftImag [i]);

                int L2 = 1 + sampleCount / 2;

                for (int i = L2; i<sampleCount; i++)
                    pts.Add (new Point (frequencyScale [i], magnitude [i]));

                for (int i = 0; i<L2; i++)
                    pts.Add (new Point (frequencyScale [i], magnitude [i]));

                LineView lv = new LineView (pts);
                figure.Plot (lv);
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception in PlotSpectrum: " + ex.Message);
            }
        }

        //**************************************************************************

        private void PlotSignal ()
        {
            try
            { 
                List<Point> pts1 = new List<Point> ();

                for (int i = 0; i<inputSignal.Length; i++)
                    pts1.Add (new Point (sampleTime [i], inputSignal [i]));

                LineView lv1 = new LineView (pts1);
                lv1.Color = Brushes.Red;
                figure.Plot (lv1);
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception: " + ex.Message);
            }
        }
    }
}
