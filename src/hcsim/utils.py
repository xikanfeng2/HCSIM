import os
import sys
import random
import datetime
import subprocess as sp
import numbers

# check part
def check_exist(**params):
    """Check that files are exist as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not os.path.exists(params[p]):
            raise ValueError(
                "{} file or directory {} does not exist.".format(p, params[p]))

def check_positive(**params):
    """Check that parameters are positive as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] <= 0:
            raise ValueError(
                "Expected {} > 0, got {}".format(p, params[p]))

def check_lt_zero(**params):
    """Check that parameters are larger than zero as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] < 0:
            raise ValueError(
                "Expected {} >= 0, got {}".format(p, params[p]))


def check_int(**params):
    """Check that parameters are integers as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not isinstance(params[p], numbers.Integral):
            raise ValueError(
                "Expected {} integer, got {}".format(p, params[p]))


def check_bool(**params):
    """Check that parameters are bools as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] is not True and params[p] is not False:
            raise ValueError(
                "Expected {} boolean, got {}".format(p, params[p]))


def check_between(v_min, v_max, **params):
    """Checks parameters are in a specified range

    Parameters
    ----------

    v_min : float, minimum allowed value (inclusive)

    v_max : float, maximum allowed value (inclusive)

    params : object
        Named arguments, parameters to be checked

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] < v_min or params[p] > v_max:
            raise ValueError("Expected {} between {} and {}, "
                             "got {}".format(p, v_min, v_max, params[p]))

def check_in(choices, **params):
    """Checks parameters are in a list of allowed parameters
    Parameters
    ----------
    choices : array-like, accepted values
    params : object
        Named arguments, parameters to be checked
    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] not in choices:
            raise ValueError(
                "{} value {} not recognized. Choose from {}".format(
                    p, params[p], choices))

def randomSNPList(chrom_sizes, snp_ratio):
    snpList = {}
    for chrom, chrom_len in chrom_sizes.items():
        snpList[chrom] = {snp : random.sample(['A','T','C','G'], 2) for snp in random.sample(range(1, chrom_len+1), int(round(chrom_len * snp_ratio)))}
    return snpList

def parseSNPList(snpfile):
    snpList = {}
    with open(snpfile, 'r') as input:
        for line in input:
            info = line.strip().split('\t')
            chrom = info[0]
            pos = info[1]
            allele1 = info[2].upper()
            allele2 = info[3].upper()
            if allele1 not in ['A', 'C', 'G', 'T'] or allele2 not in ['A', 'C', 'G', 'T']:
                continue
            if chrom not in snpList:
                snpList[chrom] = {}
            snpList[chrom][pos] = (allele1, allele2)
    return snpList

def parseIgnoreList(ignorefile):
    ignorelist = []
    with open(ignorefile, 'r') as input:
        for line in input:
            if line != '':
                ignorelist.append(line.strip())
    return ignorelist

def random_cnv():
    # 定义数字和对应的概率e
    numbers = [0, 2, 3, 4, 5]
    probabilities = [0.2, 0.3, 0.25, 0.15, 0.1]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_mirrored_cnv():
    # 定义数字和对应的概率e
    numbers = [2, 3, 4, 5, 6, 7, 8, 9]
    probabilities = [0.05, 0.2, 0.3, 0.2, 0.1, 0.05, 0.05, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_WGD():
    # 定义数字和对应的概率e
    numbers = [2, 3, 4, 5]
    probabilities = [0.5, 0.3, 0.15, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_CNL():
    # 定义数字和对应的概率e
    numbers = [1, 2, 3, 4, 5]
    probabilities = [0.4, 0.3, 0.15, 0.1, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def assign_cells_to_clones(cell_no, clone_no):
    # # Ensure at least one cell for each clone
    # clones = [1] * clone_no
    # remaining_cells = cell_no - clone_no

    # # Assign remaining cells to clones
    # for _ in range(remaining_cells):
    #     clone_index = random.randint(0, clone_no - 1)
    #     clones[clone_index] += 1

    # return clones
    # Initialize a list to hold the number of cells assigned to each clone
    clones = [0] * clone_no
    
    # Distribute the cells to the clones
    for i in range(cell_no):
        # Determine which clone to assign the current cell to
        clone_index = i % clone_no
        # Increment the cell count for the determined clone
        clones[clone_index] += 1
    
    return clones

def random_sampling(m, n, k):
    # m for total reads no, n for cell no, k for per reads no in cell
    if n * k < m:
        # Randomly sample without replacement
        sampled_lines = random.sample(range(m), n * k)
    else:
        # Randomly sample with replacement
        sampled_lines = [random.randint(0, m - 1) for _ in range(n * k)]

    # Reshape the sampled lines into a 2D array
    sampled_lines_2d = [sorted(sampled_lines[i:i+k]) for i in range(0, len(sampled_lines), k)]
    
    # add three lines to each line index
    new_sampled_lines_2d = []
    for lines in sampled_lines_2d:
        temp = []
        for line in lines:
            temp += [line*4+1, line*4+2, line*4+3, line*4+4]
        new_sampled_lines_2d.append(temp)
    return new_sampled_lines_2d

def root_path():
    return os.path.dirname(os.path.abspath(__file__))

class bcolors:
    HEADER = '\033[95m'  # 用于高亮显示标题或头部信息的颜色代码，通常为浅紫色。
    OKBLUE = '\033[94m'  # 一种蓝色，可能用于表示正常或信息性的消息。
    BBLUE = '\033[96m'  # 另一种蓝色，可能用于表示警告或其他需要引起注意的信息。
    OKGREEN = '\033[92m'  # 用于表示成功或正面的消息，通常为绿色。
    WARNING = '\033[93m'  # 用于表示警告或需要注意的信息，通常为黄色。
    FAIL = '\033[91m'  # 用于表示错误或失败的信息，通常为红色。
    ENDC = '\033[0m'  # 用于重置终端文本颜色，使其恢复到默认颜色。
    BOLD = '\033[1m'  # 用于加粗文本。
    UNDERLINE = '\033[4m'  # 用于下划线


# ProgressBar 在命令行界面显示进度条，向用户提供任务完成的百分比和状态更新。
class ProgressBar:

    def __init__(self, total, length, lock=None, counter=0, verbose=False, decimals=1, fill=chr(9608),
                 prefix='Progress:', suffix='Complete'):
        """
        :param total: 任务的总步骤数。
        :param length: 进度条的长度，以字符为单位。
        :param lock: 用于线程安全的锁，如果进度条在多线程环境中使用，需要这个参数。
        :param counter: 当前任务完成的步骤数，默认为0。
        :param verbose: 是否显示详细的信息，默认为False。
        :param decimals: 进度百分比的小数位数。
        :param fill: 进度条中填充部分的字符，默认为 █ ,这是一个标准的块状字符
        :param prefix: 进度条前缀，默认为"Progress:"。
        :param suffix: 进度条后缀，默认为"Complete"。
        """
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.lock = lock
        self.counter = counter
        assert lock is not None or counter == 0
        self.verbose = verbose

    def progress(self, advance=True, msg=""):
        """
        更新进度条的状态。
        :param advance: 是否增加计数器，进度增加，默认为True。
        :param msg: 要显示的附加信息。
        """
        if self.lock is None:
            self.progress_unlocked(advance, msg)
        else:
            self.progress_locked(advance, msg)
        return True

    def progress_unlocked(self, advance, msg):
        """
        更新进度条的状态，不使用锁。
        :param advance: 是否增加计数器，进度增加，默认为True。
        :param msg: 要显示的附加信息。
        :return:
        """
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            self.counter += 1
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter / float(self.total)))
        filledLength = int(self.length * self.counter // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = f'{self.prefix} |{bar}| {percent}% {self.suffix}'
        msg = f'[{datetime.datetime.now():%Y-%b-%d %H:%M:%S}] {msg}'
        if not self.verbose:
            toprint = rewind + result + f" [{msg}]"
        else:
            toprint = rewind + msg + "\n" + result
        write(toprint)
        flush()
        if self.counter == self.total:
            write("\n")
            flush()

    def progress_locked(self, advance, msg):
        """
        更新进度条的状态，*使用锁。
        :param advance: 是否增加计数器，进度增加，默认为True。
        :param msg: 要显示的附加信息。
        :return:
        """
        flush = sys.stderr.flush
        write = sys.stderr.write
        # 1. 更新操作
        if advance:
            # 1.1.使用 self.counter.get_lock() 获取锁，确保线程安全地更新计数器。
            with self.counter.get_lock():
                # 1.2.进度条的当前进度
                self.counter.value += 1
        # 2.计算进度条的完成百分比，使用格式化字符串将结果保留到指定的小数位数。
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter.value / float(self.total)))
        # 3.计算进度条中填充部分的长度，使用整数除法 // 将计数器的值除以总步骤数，然后乘以进度条的长度。
        filledLength = int(self.length * self.counter.value // self.total)
        # 4.bar 根据 filledLength 构建进度条的字符串表示：
        # 使用 self.fill 指定的字符填充已完成部分，使用 '-' 表示未完成部分。
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        # # 5.rewind 是一个 ANSI 转义序列，用于清除当前行并返回行首，这样每次更新进度条时都能在同一行显示
        rewind = '\x1b[2K\r'
        # 6.result 构建进度条的字符串，包括前缀、进度条本身、完成百分比和后缀。
        # msg 格式化消息字符串，包括当前的时间戳和传入的 msg 参数。
        result = f'{self.prefix} |{bar}| {percent}% {self.suffix}'
        msg = f'[{datetime.datetime.now():%Y-%b-%d %H:%M:%S}] {msg}'
        # 7.根据 self.verbose 的值，决定是否在进度条前显示额外的详细信息。
        # 如果False，只显示进度条和进度条提示消息；如果True，先显示工作函数的处理消息，然后是进度条。
        if not self.verbose:
            toprint = rewind + result + f" [{msg}]"
        else:
            toprint = rewind + msg + "\n" + result
        # 8.使用 self.lock 确保在写入和刷新时不会与其他线程冲突。
        with self.lock:
            sys.stderr.write(toprint)
            sys.stderr.flush()
            # 9.如果 self.counter.value 等于 self.total，表示任务已经完成。
            # 在任务完成后，输出一个换行符并刷新标准错误流，以确保进度条后面没有残留的内容。
            if self.counter.value == self.total:
                write("\n")
                flush()

# def error(msg):
#     log(msg=msg, level="ERROR")
#     sys.exit(0)

def which(program):
    """
    查找并返回一个可执行文件的路径。如果找到了，它就返回该文件的路径；如果没有找到，它就返回None。
    :param program:
    :return:
    """
    import os
    def is_exe(fpath):
        # 检查给定的文件路径是否存在并且具有执行权限。
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    # 提供了完整路径（fpath不为空）
    if fpath:
        # 该路径下的文件是否存在且可执行
        if is_exe(program):
            return program
    else:
        # 函数会搜索环境变量PATH中列出的所有目录
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            # 寻找第一个匹配且可执行的文件。
            if is_exe(exe_file):
                return exe_file

    return None

def runcmd(cmd, log):
    '''
    simple version of run command
    run shell command and write error log
    '''
    with open(log, 'a') as errlog:
        try:
            # Run the command
            proc = sp.Popen(
                cmd,
                shell=True,  # Allows the use of shell-specific syntax
                stdout=sp.PIPE,  # Captures standard output
                stderr=errlog  # Redirects standard error to the log file
            )
            proc.wait()  # Wait for the process to complete
            
        except Exception as e:
            print(f"Error executing command: {cmd}\n{e}")

# backup
# def runcmd(cmd, xdir, out=None, log="log"):
#     """
#     在指定目录下执行一个命令，并将输出和错误日志保存到文件中。
#     :param cmd: 要执行的命令字符串。
#     :param xdir: 执行命令的目录路径。
#     :param out: 选参数，指定输出文件的名称。如果为 None，则不保存标准输出。
#     :param log: 日志文件的名称，默认为 "log"。
#     :return:
#     """
#     j = os.path.join
#     tmp = log + '_TMP'


#     # 1. 确定标准输出文件，如果 out 参数为 None，则使用 os.devnull
#     sout_path = j(xdir, out) if out else os.devnull
#     tmp_path = j(xdir, tmp)
#     log_path = j(xdir, log)

#     # 2. 使用 with open 自动管理文件句柄
#     with open(sout_path, 'a') as sout, open(tmp_path, 'a') as serr:
#         # 使用 subprocess.Popen 启动子进程执行命令
#         proc = sp.Popen(cmd, shell=True, stdout=sout, stderr=sp.PIPE)

#         # 逐行读取标准错误输出并写入到 stderr 和临时文件
#         for line in iter(proc.stderr.readline, b''):  # readline 返回字节数据
#             sys.stderr.write(line.decode())  # 将字节解码为字符串并写入标准错误
#             serr.write(line.decode())  # 将标准错误写入临时文件

#         # 等待子进程完成
#         proc.stderr.close()
#         proc.wait()

#     # 读取临时文件内容并清理 ANSI 转义序列
#     with open(tmp_path, 'r') as i, open(log_path, 'a') as o:
#         for line in i:
#             if 'Progress' not in line:
#                 # 移除 ANSI 转义序列
#                 o.write(re.sub(r'\033\[[0-9]*m', '', line))

#     # 删除临时文件
#     os.remove(tmp_path)
