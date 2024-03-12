using System;
using System.Collections.Generic;
using System.Numerics;
using System.Windows;
using System.Windows.Media;

using MathNet.Numerics.IntegralTransforms;
using Plot2D_Embedded;

namespace ComplexFFT
{
    public partial class MainWindow : System.Windows.Window
    {
        const int sampleCount = 1024;
        const double sampleFreq = 2048; // samples per second
        const double resolution = sampleFreq / sampleCount;

        Complex [] inputSignal = new Complex [sampleCount];
        Complex [] workBuffer  = new Complex [sampleCount]; // FFTs are in-place
        double  [] sampleTime  = new double  [sampleCount];

        List<double> frequencies = new List<double> () { 4 * resolution, 16 * resolution };//, 200, 400 };
   //   List<double> frequencies = new List<double> () {23.456, 65.432 };
        List<double> amplitudes = new List<double> ()  { 4, 1, 0.5 }; // one for each frequency

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
            //
            // sample times
            //
            for (int i = 0; i<sampleCount; i++)
                sampleTime [i] = i / sampleFreq;

            //
            // complex input signal
            //
            for (int i=0; i<sampleCount; i++)
            {
                double re = 0, im = 0;

                for (int j=0; j<frequencies.Count; j++)
                {
                    re += amplitudes [j] * Math.Cos (2 * Math.PI * frequencies [j] * sampleTime [i]);
                    im += amplitudes [j] * Math.Sin (2 * Math.PI * frequencies [j] * sampleTime [i]);
                }

                inputSignal [i] = new Complex (re, im);
            }

            //
            // Copy input signal to work buffer. FFT is in-place
            //
            inputSignal.CopyTo (workBuffer, 0);

            Fourier.Forward (workBuffer);
        }

        //*****************************************************************************

        private void PlotSelectCB_Click (object sender, RoutedEventArgs e)
        {
            figure.Clear ();

            if (PlotSignalCB.IsChecked == true) PlotSignal ();
            if (PlotSpectrumCB.IsChecked == true) PlotSpectrum ();
            //if (PlotFilterCB.IsChecked == true) PlotFilterCoefficients ();

            figure.RectangularGridOn = true;
        }

        //*****************************************************************************

        // swap halves

        private void PlotSpectrum ()
        {
            List<Point> pts = new List<Point> ();

            double [] frequencyScale = Fourier.FrequencyScale (workBuffer.Length, sampleFreq);
            double [] magnitude = new double [workBuffer.Length];

            for (int i=0; i<workBuffer.Length; i++)
                magnitude [i] = Math.Sqrt (workBuffer [i].Real * workBuffer [i].Real + workBuffer [i].Imaginary * workBuffer [i].Imaginary);

            for (int i=0; i<workBuffer.Length; i++)
                magnitude [i] = 20 * Math.Log10 (magnitude [i]);

            int L2 = 1 + workBuffer.Length / 2;

            for (int i = L2; i<workBuffer.Length; i++)
                pts.Add (new Point (frequencyScale [i], magnitude [i]));

            for (int i = 0; i<L2; i++)
                pts.Add (new Point (frequencyScale [i], magnitude [i]));

            LineView lv = new LineView (pts);
            figure.Plot (lv);
        }

        //**************************************************************************

        private void PlotSignal ()
        {
            List<Point> pts1 = new List<Point> ();
            List<Point> pts2 = new List<Point> ();

            for (int i = 0; i<inputSignal.Length; i++)
            {
                pts1.Add (new Point (sampleTime [i], inputSignal [i].Real));
                pts2.Add (new Point (sampleTime [i], inputSignal [i].Imaginary));
            }

            LineView lv1 = new LineView (pts1);
            LineView lv2 = new LineView (pts2);
            lv1.Color = Brushes.Red;
            lv2.Color = Brushes.Green;
            figure.Plot (lv1);
            figure.Plot (lv2);
        }
    }
}

