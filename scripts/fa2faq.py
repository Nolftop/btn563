#fastq generator coz I dont have RAW reads
import random
import string

random.seed(563)

def id2header(id):
    header = '@NIKITA:{}:HW2:{}:{}:{}:{} {}:N:18:{}'.format(id.split('|')[1],
                                                                     str(random.randint(1, 9)),
                                                                     str(random.randint(1000, 9999)),
                                                                     str(random.randint(10000, 99999)),
                                                                     str(random.randint(100000, 999999)),
                                                                     str(random.randint(1, 9)),
                                                                     str(random.randint(1, 9)))
    return header

def seq2q(nuc):
    quality = '890ABCDEFI=<>'
    return ''.join(random.choice(quality) for i in range(len(nuc)))


if __name__ == '__main__':
    with open('../data/ncbi_query/mind.fasta', 'r') as fasta, open('../data/fastq/mind.fastq', 'w') as fastq:
        for line in fasta:
            if line.startswith('>'):
                header = id2header(line)
                one_line_dna = ''
            if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
                one_line_dna += line.strip('\n')
            if line.startswith('\n'):
                fastq.write(header)
                fastq.write('\n')
                fastq.write(one_line_dna)
                fastq.write('\n')
                fastq.write('+')
                fastq.write('\n')
                fastq.write(seq2q(one_line_dna))
                fastq.write('\n')
                fastq.write('\n')
