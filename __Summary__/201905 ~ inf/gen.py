# python -u gen.py
import os
def getPath():
    return os.path.dirname(os.path.abspath(__file__))

def checkSuffix(fileName):
    return os.path.splitext(fileName)[-1][1:] == 'md'

def getFiles():
    for root, dirs, files in os.walk(getPath()):
        return list(filter(checkSuffix, [os.path.join(root, f) for f in files]))

def main():
    os.chdir(getPath()) # change path to path of this python source
    files = sorted(getFiles()) # sort files by alphabetical order
    newFile = ''
    for file in files:
        with open(file, 'r', encoding='utf-8') as f:
            newFile += f.read()
        newFile += '\n\n'
    if not os.path.exists('generate'):
        os.makedirs('generate')
    with open('generate/summary.md', 'w', encoding='utf-8') as f:
        f.write(newFile)

if __name__ == '__main__':
    main()
