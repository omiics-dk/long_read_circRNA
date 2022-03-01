import argparse
import sys
import os
from pathlib import Path
import subprocess

class ArgumentCaller():
    
    __version__ = 2
    
    def __init__(self):
        print(
"""\033[1m  
  _                    ___             _      _        ___ _  _   _   
 | |   ___ _ _  __ _  | _ \___ __ _ __| |  __(_)_ _ __| _ \ \| | /_\  
 | |__/ _ \ ' \/ _` | |   / -_) _` / _` | / _| | '_/ _|   / .` |/ _ \ 
 |____\___/_||_\__, | |_|_\___\__,_\__,_| \__|_|_| \__|_|_\_|\_/_/ \_\ 
               |___/                                                    

Version: v{}\033[0m""".format(self.__version__))
        
        parser = argparse.ArgumentParser(
            description="Description: CircRNA detection in nanopore data",
            usage='''long_read_circRNA <command> [args]
            
Available subcommands:
    run                 Main command for finding circRNAs in long read nanopore data
    download-data       Download any required data for references and prepare
    check-installation  Check if all of the required software is available
''')
        parser.add_argument('command', help='Specify which subcommand should be used')
        args = parser.parse_args(sys.argv[1:2])
        command = args.command
        # Correct for inputs with - in the name
        command = command.replace("-","_")
        if not hasattr(self, command):
            print("Unregcognized command")
            parser.print_help()
            sys.exit(1)
        getattr(self, command)()
        
    def run(self):
        parser = argparse.ArgumentParser(
            description="Run the program for detecting circRNAs in nanopore data",
            usage="long_read_circRNA run sample [args]"
        )
        
        parser.add_argument('sample', help="Provide a sample input .fq.gz file that should be processed by the tool")
        parser.add_argument('--reference-path', default="./data", help="Provide a path for where the reference data is located.")
        parser.add_argument('--species', default='human', choices=['human', 'mouse'])
        
        args = parser.parse_args(sys.argv[2:])
    
        reference_path = args.reference_path
        species = args.species
        
        # Check if reference_path exists
        if not os.path.exists(reference_path):
            raise Exception("'{}' does not exists! Please make sure that the path is written correctly".format(reference_path))
        # Check if reference_path is a directory:
        if not os.path.isdir(reference_path):
            raise Exception("'{}' is not a directory! Please provide a directory that contains the reference data".format(reference_path))
        # Check if reference_path contains species in the end
        if os.path.basename(reference_path)==species:
            raise Exception("Looks like you have provided the direct path to the species reference directory '{}'. Please simply provide data directory '{}'".format(
                reference_path, os.path.dirname(reference_path)
            ))
        # Check if reference_path + species exists
        if not os.path.exists(os.path.join(reference_path, species)):
            raise Exception("The species subpath {} in the directory {} does not exists".format(species, reference_path))
        # Check if the reference_path + species is a directory
        if not os.path.isdir(os.path.join(reference_path, species)):
            raise Exception("The species subpath {} in the directory {} is not a directory!".format(species, reference_path))
        
        # Check if sample exists
        if not os.path.exists(args.sample):
            raise Exception("Sample file '{}' does not exist!".format(args.sample))
        # Check if sample ends with ".fq.gz"
        if not args.sample.endswith(".fq.gz"):
            raise Exception("Sample file '{}' does not end with '.fq.gz'".format(args.sample))
        # Prepare sample_path and sample_name
        reference_path = str(Path(reference_path).resolve())
        sample_path = str(Path(os.path.dirname(args.sample)).resolve())
        sample_name = os.path.basename(args.sample).replace('.fq.gz', '')
        
        if not os.path.exists(os.path.join(os.getcwd(), "scripts")) and not os.path.exists(os.path.join(os.getcwd(), "long_read_circRNA.py")):
            raise Exception("Please only run this script in the original repository .")
        
        print()
        print("\033[1m Starting process with sample: {}\033[0m".format(sample_name))
        print("\033[1m Reference path: {}\033[0m".format(reference_path))
        print("\033[1m Sample path: {}\033[0m".format(sample_path))
        print("\033[1m Species: {}\033[0m".format(sample_name))
        
        original_directory = os.getcwd()
        
        # Main process for circRNA detection
        subprocess.run(["bash", "scripts/blat_nanopore_v5.5.sh", sample_path, sample_name, species, reference_path])
        
        print("\033[1mcircRNA detection has finished\033[0m")
        print("\033[1mStarting the novel exon and alternative usage script\033[0m")
        
        os.chdir(original_directory)
        
        subprocess.run(["bash", "scripts/novel_exons_and_alternative_usage_v7.0.sh", sample_name, species, reference_path])
    
    def download_data(self):
        parser = argparse.ArgumentParser(
            description="Download some test data and required reference data for human and mouse",
            usage="long_read_circRNA download-data --data-output ./data --test-output ./test_fastq [args]"
        )
        
        parser.add_argument("--data-output", default="./data", help="Where the reference data should be stored")
        parser.add_argument("--test-output", default="./test_fastq", help="Location for the test data")
        parser.add_argument("--skip-human", action="store_true", help="Skip downloading the human reference files")
        parser.add_argument("--skip-mouse", action="store_true", help="Skip downloading the mouse reference files")
        parser.add_argument("--skip-testdata", action="store_true", help="Skip downloading the test data")
        args = parser.parse_args(sys.argv[2:])
        
        skip_check = {True: "skip", False: "no"}
        
        skip_human = skip_check[args.skip_human]
        skip_mouse = skip_check[args.skip_mouse]
        skip_testdata = skip_check[args.skip_testdata]
        
        call = subprocess.run(["bash", "get_required_data.sh", args.data_output, args.test_output,
                               skip_testdata, skip_human, skip_mouse])
    
    def check_installation(self):
        parser = argparse.ArgumentParser(
            description="Check if environment for running the script is properly setup by checking if all of the required software is available",
            usage="long_read_circRNA check-installation [args]"
        )
        
        args = parser.parse_args(sys.argv[2:])
        
        tools = {"bedtools": "bedtools", 
                 "NanoFilt": "NanoFilt -h", 
                 "pblat": "pblat", 
                 "perl": "perl -v",
                 "samtools": "samtools --help"}
        
        failed = False
        
        print()
        
        for tool in tools:
            print("\033[1mChecking for {}\033[0m".format(tool))
            try:
                output = subprocess.run(tools[tool].split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5)
                if output.returncode!=0:
                    print("\033[1;93m\tThere might be a problem with {}\033[0m".format(tool))
                    failed = True
            except FileNotFoundError:
                print("\033[1;91m\tUnable to find {}!\033[0m".format(tool))
                failed = True
            except Exception as e:
                print(e)
        
        if failed:
            print()
            print("\033[1;91mERROR: Some of the requirement software is either missing or has potential problems!\033[0m")
            print()
            sys.exit(1)
        else:
            print()
            print("\033[1mAll of the expected software requirements are present!\033[0m")
            print()        
        
if __name__=="__main__":
    ArgumentCaller()