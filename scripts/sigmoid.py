import sys
import math

if __name__ == "__main__":

        x = float(sys.argv[1]) #should be a score from 0.0 to 1. #should be a score from 0.0 to 1.00
        output_filename = sys.argv[2]

        y = 0.49 / (0.5 + math.exp(-2.0*(x-0.3)))

        output_file = open(output_filename, 'w')
        output_file.write("%f\n" % y)
        output_file.close()

