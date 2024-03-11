using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Media;

using MathNet.Filtering.FIR;

using Plot2D_Embedded;
//using System.Runtime.Remoting.Metadata.W3cXsd2001;

namespace FirFilterPlots
{
    public partial class MainWindow : Window
    {
        double [] filterCoefs;
        private OnlineFirFilter filter;

        readonly double sampleRate = 4096; 
        readonly double cutoff = 512;
        readonly int decimation = 2;
        readonly List<double> frequencies = new List<double> () {512};
        readonly double duration = 0.1; // seconds
        readonly double plotDuration = 0.1;

        readonly List<double> inputSampleTimes = new List<double> ();
        readonly List<double> inputSignal = new List<double> ();

        double [] outputSignal;
        double [] outputSampleTimes;

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

        //****************************************************************************************
        //
        // make input signal, filter coefficients and output signal
        //
        private void DoCalculations ()
        {
            try
            {
                // input signal
                for (double t = 0; t<duration; t+=1/sampleRate)
                    inputSampleTimes.Add (t);

                foreach (double t in inputSampleTimes)
                {
                    double s = 0;
                    foreach (double f in frequencies)
                        s += Math.Sin (2 * Math.PI * f * t);
                    inputSignal.Add (s);
                }

                // generate filter coefficients
                filterCoefs = MathNet.Filtering.FIR.FirCoefficients.LowPass (sampleRate, cutoff); 
                filter = new OnlineFirFilter (filterCoefs);
                Console.WriteLine (filterCoefs.Length + " filter coefficients");

                // run filter to generate output
                double [] fullOutputSignal = filter.ProcessSamples (inputSignal.ToArray ());

                // decimate
                outputSignal      = new double [(int) Math.Ceiling ((double) fullOutputSignal.Length / decimation)];
                outputSampleTimes = new double [(int) Math.Ceiling ((double) fullOutputSignal.Length / decimation)];

                for (int i = 0; i<outputSignal.Length; i++)
                {
                    outputSignal      [i] = fullOutputSignal [i * decimation];
                    outputSampleTimes [i] = inputSampleTimes [i * decimation];
                }
            }

            catch (Exception ex)
            {
                Console.WriteLine ("Exception: " + ex.Message);
            }
        }

        //************************************************************************
        //
        // Plot input signal
        //
        private void PlotInputSignal ()
        {
            List<Point> pts = new List<Point> ();

            for (int i = 0; i<inputSampleTimes.Count; i++)
                pts.Add (new Point (inputSampleTimes [i], inputSignal [i]));

            LineView lv = new LineView (pts);
            figure.Plot (lv);

            // limit display for better view of first part
            figure.GetAxes (out _, out _, out double ymin, out double ymax);
            figure.SetAxes (0, plotDuration, ymin, ymax);

            figure.RectangularGridOn = true;
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

            figure.Plot (lv);
            figure.RectangularGridOn = true;
        }

        //***********************************************************************
        //
        // Plot filter output
        //
        private void PlotOutputSignal ()
        {
            List<Point> pts = new List<Point> ();

            for (int i = 0; i<outputSampleTimes.Length; i++)
                pts.Add (new Point (outputSampleTimes [i], outputSignal [i]));

            LineView lv = new LineView (pts);
            lv.Color = Brushes.Red;
            figure.Plot (lv);

            // limit display for better view of first part
            figure.GetAxes (out _, out _, out double ymin, out double ymax);
            figure.SetAxes (0, plotDuration, ymin, ymax);

            figure.RectangularGridOn = true;
        }

        //*****************************************************************************

        private void PlotSelectCB_Click (object sender, RoutedEventArgs e)
        {
            figure.Clear ();

            if (PlotInputCB.IsChecked == true) PlotInputSignal ();
            if (PlotOutputCB.IsChecked == true) PlotOutputSignal ();
            if (PlotFilterCB.IsChecked == true) PlotFilterCoefficients ();

            figure.RectangularGridOn = true;
        }
    }
}
