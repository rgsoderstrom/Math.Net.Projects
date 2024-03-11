using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//using PerformanceAnalysis.FFT;
//using PerformanceAnalysis.LinearAlgebra;
//using MathNet.Numerics.Transformations;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Filtering;
using System.IO;

namespace FirFilter
{
    internal class MainProgram
    {
        static void Main (string [] args)
        {
            double [] coefs = MathNet.Filtering.FIR.FirCoefficients.LowPass (19150, 7500);

            Console.WriteLine (coefs.Length);

            //    foreach (double d in coefs)
            //       Console.WriteLine (d);


            MathNet.Filtering.FIR.OnlineFirFilter filter = new MathNet.Filtering.FIR.OnlineFirFilter (coefs);

            double [] input = new double [100];

            double[] output = filter.ProcessSamples (input);

            PlotWindow win = new PlotWindow ();

            Console.ReadKey ();
        }
    }
}
