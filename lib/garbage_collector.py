import os,glob,sys

class GarbageCollector:

    file_list = None
    dir_list = None

    def __init__(self):
        file_list = []
        dir_list = []
        print('creating garbage collector')

    # TODO: adds to list of files to remove
    def add_file_to_garbage_list(self,file_path):
        print('adding to file list')
        if not os.path.exists(file_path):
            sys.stderr.write(f"Error, {file_path} does not exist, not adding to file list\n")
            return False
        else:
            self.file_list.append(file_path)
            return True


    # TODO: adds to list of directories to remove
    def add_directory_to_garbage_list(self,dir_path):
        print('adding to dir list')
        if not os.path.exists(dir_path):
            sys.stderr.write(f"Error, {dir_path} does not exist, not adding to dir list\n")
            return False
        else:
            self.dir_list.append(dir_path)
            return True

    def remove_files(self):
        for file_path in self.file_list:
            if not os.path.exists(file_path):
                sys.stderr.write(f"Error, {file_path} does not exist\n")
                continue
            try:
                os.remove(file_path)
            except Exception as e:
                sys.stderr.write(f"Error removing {file_path}:\n{e}\n")

    # TODO: should I remove non-empty directories if they are added?
    def remove_directories(self):
        for dir_path in self.dir_list:
            if not os.path.exists(dir_path):
                sys.stderr.write(f"Error, {dir_path} does not exist\n")
                continue
            if not os.listdir(dir_path):
                sys.stdout.write(f"{dir_path} is empty, removing\n") 
                try:
                    os.rmdir(dir_path)
                except Exception as e:
                    sys.stderr.write(f"Error removing {dir_path}:\n{e}\n")
            else:
                sys.stdout.write(f"{dir_path} is not empty, skipping\n")
                continue
