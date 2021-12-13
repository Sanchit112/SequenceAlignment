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

    def align2(self, x, y):
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

        return ''.join(xans[i:l + 1]), ''.join(yans[i:l + 1]), dp[m][n]

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
        x, y, cost = self.align(X, Y)
        m1x, m1y = tracemalloc.get_traced_memory()
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
