# Long read circRNA

This is a tool designed to be run on a linux system since it uses many standard GNU tools.
It might also work on Mac OS X system.

Provides a workflow for detection of circRNA in nanopore data.

## Installation

For conveince clone the repository to the home directory:

```
cd ~
git clone https://github.com/omiics-dk/long_read_circRNA.git
```

Alternatively you can clone the repository somewhere else and make a symbolic link to the repository.

```
ln -s path/to/long_read_circRNA ~/long_read_circRNA
```

If you have [conda](https://docs.conda.io/en/latest/) you can setup an environment with all of the required dependencies.

```
conda env create -f long_read_circRNA.yml
```

Or you can build an environment by your self.

```
conda create -n long_read_circRNA -c bioconda bedtools=2.29.2 samtools pblat nanofilt deeptools
```

If perl is not installed on your system you can include perl in the installation:

```
conda create -n long_read_circRNA -c bioconda -c conda-forge perl bedtools=2.29.2 samtools pblat nanofilt deeptools
```

To use the environment you need to activate it afterwards.

```
conda activate long_read_circRNA
```

To make the long_read_circRNA script available in the commandline everywhere in the system you can place the file somewhere accessible in the PATH, such as /usr/local/bin or ~/bin

## Command line interface
 
Working with the tools is done through using the `long_read_circRNA` python CLI tool. The tool includes three different subcommands can be called:

* run
* download-data
* check-installation

To view the help message you can run the following command `./long_read_circRNA -h`.


```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1
usage: long_read_circRNA <command> [args]
            
Available subcommands:
    run                 Main command for finding circRNAs in long read nanopore data
    download-data       Download any required data for references and prepare
    check-installation  Check if all of the required software is available

Description: CircRNA detection in nanopore data

positional arguments:
  command     Specify which subcommand should be used

optional arguments:
  -h, --help  show this help message and exit

```

## Check installation

Before running the script, it is a good idea to use the check-installation command to see if the required programs are available.

`./long_read_circRNA check-installation -h`
```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1
usage: long_read_circRNA check-installation [args]

Check if environment for running the script is properly setup by checking if
all of the required software is available

optional arguments:
  -h, --help  show this help message and exit

```

### Working example

When running `./long_read_circRNA check-installation` you will get a summary of all the tools that are check for.
If there is no message for an individual tool then it has been found and there where no problems with running it.

When the installation is properly installed and you run the check-installation subcommand you will find that returns that all of the expected software requirments are present.

```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1

Checking for bedtools
Checking for NanoFilt
Checking for pblat
Checking for perl
Checking for samtools

All of the expected software requirements are present!
```

### Failed example

Running the `./long_read_circRNA check-installation` when some of the required tools are missing then the script will return and error and list which tools are missing. If the program returns an error when called then will return a warning that there might be an issue with the program that should be checked.
```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1

Checking for bedtools
        Unable to find bedtools!
Checking for NanoFilt
        Unable to find NanoFilt!
Checking for pblat
        Unable to find pblat!
Checking for perl
Checking for samtools
        Unable to find samtools!

ERROR: Some of the requirement software is either missing or has potential problems!
```

## Download data

For conveince using the download-data subcommand provides with necessary reference files to run the analysis on human and on mice. It also provides a small test dataset to try out and run the software with.

Running `./long_read_circRNA download-data -h` shows different arguments that can be provided.

```  
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1
usage: long_read_circRNA download-data --data-output ./data --test-output ./test_fastq [args]

Download some test data and required reference data for human and mouse

optional arguments:
  -h, --help            show this help message and exit
  --data-output DATA_OUTPUT
                        Where the reference data should be stored
  --test-output TEST_OUTPUT
                        Location for the test data
  --skip-human          Skip downloading the human reference files
  --skip-mouse          Skip downloading the mouse reference files
  --skip-testdata       Skip downloading the test data
  --installation-path INSTALLATION_PATH
                        Provide a path for where the script is installed.
                        Default is ~/long_read_circRNA
```

By default the tool will automatically download references to `./data` and the test data to `./test_fastq`. These can be changed by using `--data-output` and `--test-output` flags.
Different datasets can be skipped if they should not be downloaded, by using the skip flags.

## Run

For running the software you use the run subcommand, and provide a sample `.fq.gz` file. There are serval options for changing the behavoiur of the run command such as:

* Specifying which species you are using, currently there are only options for human and mouse. Default is human.
* Defining a path to where the reference files are location. Default is ./data.
  * It should be noted that it automatically searchs for ./data/human if you choose human as the species, and therefore you should only provide the ./data directory
* Providing the path to where the scripts are located. Defaults to the home directory `~/long_read_circRNA/scripts`
* Providing where the outputs will be written. Default is the current directory it is running it, where it will automatically create directory with the sample name

```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2.1
usage: long_read_circRNA [-h] [--reference-path REFERENCE_PATH]
                         [--species {human,mouse}] [--script-path SCRIPT_PATH]
                         [--output-path OUTPUT_PATH] [--dry-run]
                         sample

Run the program for detecting circRNAs in nanopore data

positional arguments:
  sample                Provide a sample input .fq.gz file that should be
                        processed by the tool

optional arguments:
  -h, --help            show this help message and exit
  --reference-path REFERENCE_PATH
                        Provide a path for where the reference data is
                        located. Default is './data'.
  --species {human,mouse}
                        Select which species the sample that is from, and
                        specify which species reference files should be used.
                        Default is set to human.
  --script-path SCRIPT_PATH
                        Specify where the long_read_circRNA scripts are
                        located. By default it assumes that they are located
                        at '~/long_read_circRNA/scripts'.
  --output-path OUTPUT_PATH, -o OUTPUT_PATH
                        Provide a path for where the output should be saved.
                        Default is the current directory. This will use the
                        output-path and create a directory in that path based
                        on the sample name provided.
  --dry-run             Perform all of the input checks without starting the
                        detection scripts
```

The script has numerous checks and error messages if it finds any errors in the provided inputs so the scripts do not start without everything being in place.

To start the circRNA quantification on a sample you can simply write:
``` bash
./long_read_circRNA run path_to_sample/sample_name.fq.gz
```

If you need to change any of the settings in the run then you can use the flags presented above.

When the repository is installed or available from the home directory, and the long_read_circRNA script is avaible in the PATH then a simple workflow for starting sample is:

``` bash
cd path_to_where_analysis_should_be_run

long_read_circRNA download-data

long_read_circRNA run sample_name.fq.gz
```

To run multiple samples you can easily use a for loop to run the script on multiple samples:

``` bash
cd path_to_where_analysis_should_be_run

long_read_circRNA download-data

samples = "sample_data/sample1.fq.gz sample_data/sample2.fq.gz sample_data/sample3.fq.gz"

for sample in $samples
do
  long_read_circRNA run sample_name.fq.gz
done
```

The scripts output a number of files, the primary one being:

*[sample].circRNA_candidates.annotated.txt*

This file shows all circRNAs detected in the nanopore data, listed with the highest expressed at the top.

