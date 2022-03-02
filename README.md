# Long read circRNA

This is a tool designed to be run on a linux system since it uses many standard GNU tools.
It might also work on Mac OS X system.

## Installation

Clone the repository with all of the scripts in the script directory.

If you have [conda](https://docs.conda.io/en/latest/) you can setup an environment with all of the required dependencies.

```
conda env create -f long_read_circRNA.yml
```

Or you can build an environment by your self.

```
conda create -n long_read_circRNA -c bioconda bedtools samtools pblat nanofilt
```

If perl is not installed on your system you can include perl in the installation:

```
conda create -n long_read_circRNA -c bioconda -c conda-forge perl bedtools samtools pblat nanofilt
```

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

Version: v2
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

Version: v2
usage: long_read_circRNA check-installation [args]

Check if environment for running the script is properly setup by checking if
all of the required software is available

optional arguments:
  -h, --help  show this help message and exit

```

### Working example

`./long_read_circRNA check-installation`
```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2

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

Version: v2

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

Version: v2
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
```

By default the tool will automatically download references to `./data` and the test data to `./test_fastq`. These can be changed by using `--data-output` and `--test-output` flags.
Different datasets can be skipped if they should not be downloaded, by using the skip flags.

## Run

To 

```
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v2
usage: long_read_circRNA run sample [args]

Run the program for detecting circRNAs in nanopore data

positional arguments:
  sample                Provide a sample input .fq.gz file that should be processed by the tool

optional arguments:
  -h, --help            show this help message and exit
  --reference-path REFERENCE_PATH
                        Provide a path for where the reference data is located.
  --species {human,mouse}

```

To start the circRNA quantification


The scripts output a number of files, the primary one being:

*[sample].circRNA_candidates.annotated.txt*


This file shows all circRNAs detected in the nanopore data, listed with the highest expressed at the top.

