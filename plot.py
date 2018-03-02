import numpy
from matplotlib import pyplot

def main():
    x, y = numpy.loadtxt("data.dat", unpack=True)

    pyplot.figure()
    pyplot.plot(x, y, "o")
    pyplot.show()
    pyplot.close()

if __name__ == '__main__':
    main()