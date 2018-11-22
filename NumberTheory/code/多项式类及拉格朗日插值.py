class poly(object):

    def __init__(self, *args, **kw):
        if 'args' in kw:
            args = kw['args']
        self.__n = len(args) - 1
        self.__data = []
        for i in args[::-1]:
            self.__data.append(float(i))

    def __str__(self):
        return ' + '.join('%.6f x^%d' %(j, self.__n - i) for i, j in enumerate(self.__data[::-1]))
    
    @staticmethod
    def __add(a, b):
        # a.__n > b.__n
        ret = []
        for i in range(b.__n + 1):
            ret.append(a.__data[i] + b.__data[i])
        for i in range(b.__n + 1, a.__n + 1):
            ret.append(a.__data[i])
        return poly(args=ret[::-1])

    def __add__(self, other):
        if self.__n > other.__n:
            return poly.__add(self, other)
        else:
            return poly.__add(other, self)
    
    def __mul__(self, other):
        ansn = self.__n + other.__n
        ret = [0] * (ansn + 1)
        for i1, i2 in enumerate(self.__data):
            for j1, j2 in enumerate(other.__data):
                ret[i1 + j1] += i2 * j2
        return poly(args=ret[::-1])
    
    def __sub__(self, other):
        return self.__add__(poly(args=list(map(lambda x: -x, other.__data))))

def lagrange(dic):
    ans = poly(0)
    for key, value in dic.items():
        cur = poly(value)
        for okey, ovalue in dic.items():
            if key == okey and value == ovalue:
                continue
            cur *= poly(1 / (key - okey), okey / (okey - key))
        ans += cur
    return ans

def XYtoDIC(x, y):
    ret = dict()
    for i in range(len(y)):
        ret[x[i]] = y[i]
    return ret

if __name__ == '__main__':
    x = [1, 3, 4]
    y = [2, 3, 5]
    print(lagrange(XYtoDIC(x, y)))
