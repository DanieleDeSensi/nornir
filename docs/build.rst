Building and Installing
=======================

First of all, download Nornir:

.. code-block:: shell

   $ git clone git://github.com/DanieleDeSensi/nornir.git
   $ cd nornir


To compile Nornir:

.. code-block:: shell

   $ mkdir build
   $ cd build
   $ cmake ../
   $ make

Then, you can install it:

.. code-block:: shell

   $ make install

To install it into a non-default directory *dir*, simply specify the *-DCMAKE_INSTALL_PREFIX=dir* when calling *cmake*.


Configuring the System
----------------------
Some operations performed by Nornir are hardware-specific (e.g. reading the system energy consumption, dynamically chaning the clock frequency, etc..). To perform such operations Nornir relies on the Mammut_ library, and you may want to refer to the Mammut_ documentation for more information. However, these are the main points to consider:

* To change the frequency, Mammut assumes that the *acpi_cpufreq* driver is used for scaling the clock frequency. Mammut currently does not support the *intel_pstate* driver. If your machine is using the *intel_pstate* driver, you should switch to *acpi_cpufreq*
* On Intel and AMD processors, you should load the *msr* kernel module to allow Mammut to read the energy consumption. For PowerPC and ARM processors please refer to the Mammut_ documentation.
* Changing the frequeency and reading the energy usually require *sudo* rights, and the easiest option would be to just run your application with *sudo*. If this is not possible, Mammut supports the msr-safe_ library. If you want to use msr-safe, please contact me and I will provide you the MSR registers whitelist.

.. _msr-safe: http://github.com/LLNL/msr-safe
.. _Mammut: http://danieledesensi.github.io/mammut