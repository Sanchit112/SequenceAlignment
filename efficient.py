import math
import random
import sys

from Bio import pairwise2
import tracemalloc
import time


class Solution:
    def __init__(self, alpha, delta):
        self.alpha = alpha
        self.delta = delta

    def align(self, X, Y):
        m = len(X)
        n = len(Y)

        dp = [[0 for x in range(n + 1)] for y in range(m + 1)]

        for i in range(m + 1):
            dp[i][0] = i * self.delta
        for j in range(n + 1):
            dp[0][j] = j * self.delta

        for j in range(1, n + 1):
            for i in range(1, m + 1):
                dp[i][j] = min(dp[i - 1][j - 1] + self.alpha.get(X[i - 1] + Y[j - 1]),
                               dp[i - 1][j] + self.delta,
                               dp[i][j - 1] + self.delta)

        i = m
        j = n

        x = ''
        y = ''

        while i != 0 and j != 0:
            if dp[i][j] == dp[i - 1][j - 1] + self.alpha.get(X[i - 1] + Y[j - 1]):
                i -= 1
                j -= 1
                x += X[i]
                y += Y[j]
            elif dp[i][j] == dp[i - 1][j] + self.delta:
                i -= 1
                x += X[i]
                y += '_'
            elif dp[i][j] == dp[i][j - 1] + self.delta:
                j -= 1
                x += '_'
                y += Y[j]

        return x[::-1], y[::-1], dp[m][n]

    def alignment(self, x, y):
        m = len(x)
        n = len(y)
        dp = []
        for i in range(n + m + 1):
            t = []
            for j in range(n + m + 1):
                t.append(0)
            dp.append(t)

        for i in range(n + m + 1):
            dp[i][0] = i * self.delta
            dp[0][i] = i * self.delta

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if x[i - 1] == y[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1]
                else:
                    dp[i][j] = min(dp[i - 1][j - 1] + self.alpha.get(x[i - 1] + y[j - 1]),
                                   dp[i - 1][j] + self.delta,
                                   dp[i][j - 1] + self.delta)

        l = n + m
        i = m
        j = n
        xpos = l
        ypos = l
        xans = []
        yans = []
        for k in range(l + 1):
            xans.append('')
            yans.append('')

        while not (i == 0 or j == 0):
            if x[i - 1] == y[j - 1]:
                xans[xpos] = x[i - 1]
                yans[ypos] = y[j - 1]
                xpos -= 1
                ypos -= 1
                i -= 1
                j -= 1
            elif dp[i - 1][j - 1] + self.alpha.get(x[i - 1] + y[j - 1]) == dp[i][j]:
                xans[xpos] = x[i - 1]
                yans[ypos] = y[j - 1]
                xpos -= 1
                ypos -= 1
                i -= 1
                j -= 1
            elif dp[i - 1][j] + self.delta == dp[i][j]:
                xans[xpos] = x[i - 1]
                yans[ypos] = '_'
                xpos -= 1
                ypos -= 1
                i -= 1
            elif dp[i][j - 1] + self.delta == dp[i][j]:
                xans[xpos] = '_'
                yans[ypos] = y[j - 1]
                xpos -= 1
                ypos -= 1
                j -= 1

        while xpos > 0:
            if i > 0:
                i -= 1
                xans[xpos] = x[i]
                xpos -= 1
            else:
                xans[xpos] = '_'
                xpos -= 1

        while ypos > 0:
            if j > 0:
                j -= 1
                yans[ypos] = y[j]
                ypos -= 1
            else:
                yans[ypos] = '_'
                ypos -= 1

        for i in range(l, -1, -1):
            if xans[i] == yans[i] == '_':
                i += 1
                break

        return ''.join(xans[i:l + 1]), ''.join(yans[i:l + 1])

    def start(self, fileName):
        strings = []
        s = ''

        with open(fileName, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) > 1:
                    strings.append(s)
                    s = line
                else:
                    n = len(s)
                    s = s[0:int(line) + 1] + s + s[int(line) + 1:n]
        return strings[1], s

    def spaceEfficient(self, x, y, flag):
        xLength = len(x)
        yLength = len(y)
        dp = [0 for j in range(yLength + 1)]

        for j in range(yLength + 1):
            dp[j] = j * self.delta

        for i in range(1, xLength + 1):
            prev = dp[0]
            dp[0] = i * self.delta

            for j in range(1, yLength + 1):
                mxy = self.alpha.get(x[i - 1] + y[j - 1]) + prev
                mx = dp[j] + self.delta
                my = dp[j - 1] + self.delta
                prev = dp[j]
                dp[j] = min(mx, my, mxy)

        if flag:
            return dp[-1]
        else:
            return dp

    def divide(self, X, Y):
        m = len(X)
        A = self.spaceEfficient(X[:m // 2], Y, False)
        x = X[m // 2:]
        x = x[::-1]
        y = Y[::-1]
        B = self.spaceEfficient(x, y, False)

        B = B[::-1]
        q = 0
        s = math.inf

        for k in range(len(A)):
            if A[k] + B[k] < s:
                q = k
                s = A[k] + B[k]

        return X[0:m // 2], Y[0:q], X[m // 2:], Y[q:]

    def end(self, X, Y):
        m = len(X)
        n = len(Y)

        if not X and Y:
            return '_' * n, Y, True
        if not Y and X:
            return X, '_' * m, True
        if not X and not Y:
            return '0', '0', True
        if m <= 2 or n <= 2:
            a, b = self.alignment(X, Y)
            return a, b, True
        else:
            return X, Y, False

    def conquer(self, X, Y):
        x1, y1, x2, y2 = self.divide(X, Y)

        x1, y1, e1 = self.end(x1, y1)
        x2, y2, e2 = self.end(x2, y2)

        if not e1 and not e2:
            a1, b1 = self.conquer(x1, y1)
            a2, b2 = self.conquer(x2, y2)
            return a1 + a2, b1 + b2
        elif not e1:
            a1, b1 = self.conquer(x1, y1)
            return a1 + x2, b1 + y2
        elif not e2:
            a2, b2 = self.conquer(x2, y2)
            return x1 + a2, y1 + b2
        return x1 + x2, y1 + y2

    def calculateCost(self, X, Y):
        ans = 0
        for i in range(len(X)):
            if X[i] != '_' and Y[i] != '_':
                ans += self.alpha.get(X[i] + Y[i])
            else:
                ans += self.delta
        return ans

    def bioCost(self, X, Y):
        bio = {}
        for k, v in self.alpha.items():
            bio[(k[0], k[1])] = -v

        align = pairwise2.align.globalds(X, Y, bio, -delta, -delta)
        return -1 * align[0].score

    def projectOutput(self, X, Y):
        start_time = time.time()
        tracemalloc.start()
        x, y = obj.conquer(X, Y)
        m1x, m1y = tracemalloc.get_traced_memory()
        cost = self.calculateCost(x, y)
        t = time.time() - start_time
        # print(m1x)
        # print(m1y)

        m = len(x)
        n = len(y)

        with open('output.txt', 'w') as f:
            if m <= 100:
                f.write(x)
                f.write('\n')
                f.write(y)
            else:
                f.write(x[0:50] + ' ' + x[m - 50:])
                f.write('\n')
                f.write(y[0:50] + ' ' + y[n - 50:])
            f.write('\n')
            f.write(str(cost))
            f.write('\n')
            f.write(str(t))
            f.write('\n')
            f.write(str(m1y))

    def test(self):
        obj = Solution(alpha, delta)

        with open('data.csv', 'w') as f1:
            for k in range(20):
                with open('input.txt', 'w') as f:
                    f.write('ACTG\n')
                    t = random.randint(5, 10)
                    for i in range(t):
                        f.write(str(random.randint(1, 4)) + '\n')
                    f.write('TACG\n')
                    for i in range(t):
                        f.write(str(random.randint(1, 4)) + '\n')

                tracemalloc.start()
                s1, s2 = obj.start('input.txt')
                # print(s1)
                # print(s2)
                start_time = time.time()
                a, b, cost = obj.align(s1, s2)
                m1x, m1y = tracemalloc.get_traced_memory()
                t = time.time() - start_time
                f1.write('test{},{},normal,{},{}\n'.format(i + 1, len(a) * len(b), t, m1y))
                print('{},{}'.format(t, m1y))
                tracemalloc.reset_peak()
                # print(a)
                # print(b)
                # print(obj.calculateCost(a, b))
                start_time = time.time()
                x, y = obj.conquer(s1, s2)
                m1x, m1y = tracemalloc.get_traced_memory()
                t = time.time() - start_time
                f1.write('test{},optimal,{},{}\n'.format(i + 1, len(x) * len(y), t, m1y))
                print('{},{}'.format(t, m1y))
                tracemalloc.stop()
                # print(x)
                # print(y)
                print(obj.calculateCost(x, y) == obj.bioCost(s1, s2) == cost)
                # print('*' * 100)


alpha = {'AA': 0, 'AC': 110, 'AG': 48, 'AT': 94,
         'CA': 110, 'CC': 0, 'CG': 118, 'CT': 48,
         'GA': 48, 'GC': 118, 'GG': 0, 'GT': 110,
         'TA': 94, 'TC': 48, 'TG': 110, 'TT': 0
         }

delta = 30

if len(sys.argv) < 2:
    print('Enter input file name as cmd line parameter')
else:
    obj = Solution(alpha, delta)
    s1, s2 = obj.start(sys.argv[1])
    obj.projectOutput(s1, s2)
    # obj.test()
# if len(sys.argv) < 2:
#     obj.projectOutput(True)
# elif sys.argv[1] == 'optimal':
#     obj.projectOutput(True)
# elif sys.argv[1] == 'dp':
#     obj.projectOutput(False)
# else:
#     obj.projectOutput(True)
# s1, s2 = obj.start('input.txt')
# tracemalloc.start()
# # x, y = obj.alignment(s1, s2)
# print(obj.align(s1, s2))
# m1x, m1y = tracemalloc.get_traced_memory()
# print(m1x)
# print(m1y)
