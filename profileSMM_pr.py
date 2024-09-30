#ввод-вывод с использованием BioPython из текстового файла
#ООП
#тесты
#выводим всю таблицу динамического программирования
#последовательности длиной >10

#var: 2.3 Профильная СММ, прямой алгоритм
# [f_M, f_I, f_D]

from Bio import SeqIO

class Probability_of_generation:
    def __init__(self, seq):

        #вероятности

        #Вероятности порождения для состояний совпадения
        self.probability_of_match = {'A': [0, 0.25, 0.55, 0.08],
                                     'C': [0, 0.08, 0.09, 0.33],
                                     'G': [0, 0.58, 0.18, 0.08],
                                     'T': [0, 0.08, 0.18, 0.08]
                                     }
        #Вероятности порождения для состояний вставки
        self.probability_of_insert = {'A': [0.25, 0.17, 0.25, 0.25],
                                      'C': [0.25, 0.5, 0.25, 0.25],
                                      'G': [0.25, 0.17, 0.25, 0.25],
                                      'T': [0.25, 0.17, 0.25, 0.25]
                                      }
        #Вероятности переходов между скрытыми состояниями
        self.probability_of_trans_betw_HS = {'M-M': [0.82, 0.55, 0.8, 0.9],
                                             'M-D': [0.09, 0.18, 0.1, 0.09],
                                             'M-I': [0.09, 0.27, 0.1, 0.1],
                                             'I-M': [0.33, 0.6, 0.33, 0.5],
                                             'I-D': [0.33, 0.2, 0.33, 0.33],
                                             'I-I': [0.33, 0.2, 0.33, 0.5],
                                             'D-M': [0, 0.33, 0.5, 0.5],
                                             'D-D': [0, 0.33, 0.25, 0.33],
                                             'D-I': [0, 0.33, 0.25, 0.5]}
        #вероятности выхода
        self.eps = [[0.0, 0.5, 0.3, 0.2], [0.1, 0.4, 0.2, 0.3], [0.0, 0.1, 0.3, 0.6]]

        #количество состояний для каждого блока (M, I, D)
        self.kol_cond = 4

        self.seq = seq
        self.seq_len = len(seq)
        self.matr_dp = []
        self.p = 0
        #инициализация
        self.matr_dp.append([[1], [0], [0], [0]])
        self.matr_dp.append([[0], [0], [0], [0]])
        self.matr_dp.append([[0], [0], [0], [0]])

    def prob_gen(self):

        #f_D_0 и f_M_0 для всех букв
        for i in range(1, self.seq_len):
            self.matr_dp[0][0].append(0)
            self.matr_dp[2][0].append(0)



        for i in range(self.kol_cond):  #для всех состояний
            for j in range(self.seq_len): #для всех букв
                if i != 0:

                    if j != 0:

                        #f_M
                        sub_f_M = self.probability_of_match[self.seq[j]][i] * \
                                  (self.matr_dp[0][i-1][j-1]*self.probability_of_trans_betw_HS['M-M'][i-1] +
                                   self.matr_dp[1][i-1][j-1]*self.probability_of_trans_betw_HS['I-M'][i-1] +
                                   self.matr_dp[2][i-1][j-1]*self.probability_of_trans_betw_HS['D-M'][i-1])
                        self.matr_dp[0][i].append(sub_f_M)

                    #f_D
                    sub_f_D = self.matr_dp[0][i-1][j]*self.probability_of_trans_betw_HS['M-D'][i-1] + \
                              self.matr_dp[1][i-1][j]*self.probability_of_trans_betw_HS['I-D'][i-1] + \
                              self.matr_dp[2][i-1][j]*self.probability_of_trans_betw_HS['D-D'][i-1]
                    self.matr_dp[2][i].append(sub_f_D)

                if j != 0:
                    # f_I
                    sub_f_I = self.probability_of_insert[self.seq[j]][i] * \
                                    (self.matr_dp[0][i][j - 1] * self.probability_of_trans_betw_HS['M-I'][i] +
                                    self.matr_dp[1][i][j - 1] * self.probability_of_trans_betw_HS['I-I'][i] +
                                    self.matr_dp[2][i][j - 1] * self.probability_of_trans_betw_HS['D-I'][i])
                    self.matr_dp[1][i].append(sub_f_I)

        print('Матрица динамического программирования:')
        #по всем состояниям внутри блока M
        for i in range(1, self.kol_cond):
            print('M_', i, ': ')
            for j in range(self.seq_len):
                print('x(',j, ') = ', self.matr_dp[0][i][j])
            print('\n')

        # по всем состояниям внутри блока I
        for i in range(0, self.kol_cond):
            print('I_', i, ': ')
            for j in range(self.seq_len):
                print('x(', j, ') = ', self.matr_dp[0][i][j])
            print('\n')

        # по всем состояниям внутри блока D
        for i in range(1, self.kol_cond):
            print('D_', i, ': ')
            for j in range(self.seq_len):
                print('x(', j, ') = ', self.matr_dp[0][i][j])
            print('\n')


        #по всем блокам
        for i in range(3):
            #по всем состояниям
            for j in range(self.kol_cond):
                self.p += self.matr_dp[i][j][self.seq_len - 1] * self.eps[i][j]
        print('Вероятность порождения данной СММ: ', self.p)



if __name__ == '__main__':
    file1 = open("seq1.fasta", "r")
    seq1 = SeqIO.read(file1, "fasta")
    print('Последовательность: ', seq1.seq)

    align = Probability_of_generation(seq1.seq)
    align.prob_gen()
