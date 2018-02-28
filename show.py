import numpy
from matplotlib import pyplot

def main():
    x1, x2 = numpy.loadtxt("data.dat", unpack=True)
    y1, bins1 = numpy.histogram(x1, bins=20)
    y2, bins2 = numpy.histogram(x2, bins=20)
    
    print(y1.shape, y1.shape)
    y1 = y1 / len(x1)
    y2 = y2 / len(x2)
    pyplot.figure()
    for i in range(len(y1)):
        pyplot.bar(bins1[i], y1[i], width=abs(bins1[0] - bins1[1]), color="blue")
        pyplot.bar(bins2[i], y2[i], width=abs(bins2[0] - bins2[1]), color="crimson")
    pyplot.show()

if __name__ == '__main__':
    main()