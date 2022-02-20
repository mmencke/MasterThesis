# MasterThesis

This is the source code for my Master Thesis on Credit Value Adjustment: Pricing Wrong Way Risk on Interest Rate Swaps. It is utilising functionality from QuantLib and works with the [forked repo](https://github.com/mmencke/QuantLib) (the 'MasterThesis'-branch).

The project is build in Xcode 9.2 using version 1.22 of QuantLib and version 1.75.0 of Boost. It is a standard C++ command-line tool with following settings

 | Setting                  | Value              |
 | ------------------------ | ------------------ |
 | Header Search Paths      | /usr/local/include |
 | Library Search Paths     | /usr/local/lib     |
 | Other Linker Flags       | -lQuantLib         |

The interest rate is assume to follow either a Hull-White or G2++ model (`HullWhite` and `G2` classes in QuantLib). The default intensity of the investor (bank) and the counterparty are assumed to follow the CIR++ model (`ExtendedCoxIngersollRoss`) with a constant recovery rate. CVA and DVA is calculated on a portfolio of vanilla swaps by Monte Carlo simulation using either Pseudo-Random (`MersenneTwister`) or Quasi-Random (`Sobol`) numbers. The results are exported to a csv-file and plotted in R. The code for the plots is in the Data-folder.
