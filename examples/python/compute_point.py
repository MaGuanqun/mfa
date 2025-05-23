import numpy as np

# points=np.array([37.9,32.7,76.4])


# points2 =np.array([36.6039,33.3391,76.2448])

# points=np.array([36,32,76])
# points2 =np.array([38,34,78])

# points=np.array([60,19,57])
# points2 =np.array([62,21,59])

def compute(points):

    points2=points+2
    size = np.array([143,255,255])
    
    
    

    # size=np.array([3599,1799])
    pos=points/size
    pos2=points2/size

    string = str(pos[0])+"-"+str(pos2[0])+"-"+str(pos[1])+"-"+str(pos2[1])+"-"+str(pos[2])+"-"+str(pos2[2])
    print(string)
    # print(pos)
    # print(pos2)
    print("===========")
    
    

points=[]
points.append([0.5034965034965035,0.10196078431372549,0.611764705882353])
size = np.array([143,255,255])

points.append([0.5034965034965035,0.10196078431372549,0.611764705882353])


# points.append(np.array([79, 65, 88.3,]))
pos=points[0]*size


print(pos)
# points.append(np.array([180,1266]))

# points.append(np.array([216,1284]))

# points.append(np.array([71,26,156]))

# points.append(np.array([57,23,160]))
# points.append(np.array([59,61,204]))
# points.append(np.array([72,26,156]))

# points.append(np.array([72.55,67.75,92.55]))

# for i in range(len(points)):
#     points[i] = points[i] + np.array([0.25,0.25,0.25])
# points.append(np.array([84,75,87]))
# points.append(np.array([72,67,92]))
# points.append(np.array([58,22,155]))
# points.append(np.array([50,50,203]))




# for point in points:
#     compute(point)