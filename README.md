# Repository.
Mouse sinoatrial node cell model codes.

# Background.  

This cardiomyocyte cell model was develop to assist collaboration with Manchester Medical School where 
the mouse is a popular animal model. It is biophysically detailed, i.e. it provides a spectrum of ion channel
species that regulate cell electrical activity. The AJP paper where the model was published is provided.
This work gave rise to the cell model codes in multiple languages suitable for use by experimentalists.
The model optimisation was based on a mutual information driven assessment of parameter-model output
sensitivty analysis.

# Sources and data description.  

Three sources of the same model are provided (see sources). The model is solved using Euler method,
higher order BDF formulae, and MATLAB's ode15s.

# Dependencies.  

Multiple versions of the model are provided. For the Euler method model (original code), there are no 
dependencies apart from the standard GNU C compiler. The user may wish to use GnuPlot to see the
model output.  

A more efficient and accurate solution is obtained using the Sundials library that provides implicit
solvers using BDF formulae. To use this version, Sundials (https://computing.llnl.gov/projects/sundials).
A suitably old version (around 2011) of Sundials should be used, but the latest versions may also work.  

To use the MATLAB version, MATLAB 2016 or later should be used. The MATLAB code can be compiled
using the compiler to improve performance to a certain extent.

# Install.  

No installation required.

# Use.  

Use as:
gcc -o msan c_filename.c -lm
./msan {potential arguments to the binary}

A few outputs I have used the model are as follows:  
* Kharche et al. Computing in Cardiology 2010;37:421−424.  
* Doris et al. (https://www.hh.um.es/Abstracts/Vol_34/34_11/34_11_1255.htm)  
* Kharche et al. (https://www.sciencedirect.com/science/article/pii/S2214854X22000024)  

# Uptake by other users.

This model has been uptaken by the Hund group (https://bme.osu.edu/hund-lab-excitable-cell-engineering),
Yaniv group (https://bme.technion.ac.il/en/team/assoc-prof-yael-yaniv/), and the Quinn group (https://quinnlaboratory.com/).

# Maintainer.

All versions of this models' code are maintained and supported by SR Kharche by means of the cardiac simulator package, Virtual Cardiac Physiology Laboratory.

# Provenance.

New releases of the self deployed PM3-SV components are handled by the lead developer and co-investigator Dr. SR Kharche (Lawson Health Institute, Canada), 
in consultation with Dr. Daniel Goldman (co-investigator) and Dr. C. W. McIntyre (named PI). Before each release any tests for the component 
or service are run to verify integration. New additions and features to the code base must pass all current tests 
and are approved by the lead developer. New additions by users will be fully tested and documented by the developer team. 
The code is commented to provide information regarding original author, testing, use, and maintainence.
The working repository also has documentation to assist uptake by new users.
After release, the developer team will provide support for the use of existing and new functionality. Each release will 
be accompanied by commits and changelog updates on their respective git/GitHub repositories. The third-party software 
and libraries are kept locked at specific versions. The PM3-SV components will receive updates when these libraries are 
updated for specific improvements.

# Acknowledements.

This project was generously funded by CANARIE Inc. Canada (April 2020 to March 2023). British Heart Foundation. Wellcome Trust. EPSRC. Exeter University internal funding. Western University internal funding.
The principal investigator was XYZ. The co-investigators were: ABC, DEC, GHI.
The employees were:
We also thank the Kidney Unit in Lawson, London Ontario, Canada, for hosting this project.

# Licence.

BSD 3-Clause License

Copyright (c) 2023, Sanjay R. Kharche, Ph.D.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

