import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import MetaPro_paths as mpp


if __name__ == "__main__":
    config_path = sys.argv[1]
    output_dir = sys.argv[2]

    path_obj = mpp.tool_path_obj(config_path, output_dir)