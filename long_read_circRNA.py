import argparse
import sys
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
            usage="long_read_circRNA run [args]"
        )
        
        args = parser.parse_args(sys.argv[2:])
        
        # Code for starting the detection
    
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