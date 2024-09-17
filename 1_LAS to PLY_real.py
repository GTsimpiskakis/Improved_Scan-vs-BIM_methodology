import laspy
import numpy as np
import open3d as o3d

def read_LAS_file(las_file):
    # Read LAS file
    inFile = laspy.file.File(las_file, mode='r')

    # Extract point cloud data
    points = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
    colors = np.vstack((inFile.red, inFile.green, inFile.blue)).transpose()
    colors = (colors / 65535) * 255
    return points, colors

def write_ply(points, colors, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = points[:, 0]
    vertex["y"] = points[:, 1]
    vertex["z"] = points[:, 2]
    vertex["red"] = colors[:, 0]
    vertex["green"] = colors[:, 1]
    vertex["blue"] = colors[:, 2]

    with open(filename, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for i in range(len(points)):
            f.write("{:8.3f} {:8.3f} {:8.3f} {:3.0f} {:3.0f} {:3.0f}\n".format(points[i][0], points[i][1], points[i][2], colors[i][0], colors[i][1], colors[i][2]))

# Specify input PCD file and output PLY file
las_file_path=r"C:/Python Output/Real data/LAS files/GTThesis RegPLUS 1.las"
ply_file_path = r"C:/Python Output/Real data/PLY files/SDK_real_dataset_final/PLY_SDK_EU_Real.ply"
# Read PCD and convert to PLY
points, colors= read_LAS_file(las_file_path)
write_ply(points, colors, ply_file_path) #Scale colors back to the 0-255 range for PLY
print("Conversion completed.")