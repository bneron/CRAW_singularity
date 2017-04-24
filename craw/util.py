import sys


def progress(count, total, status='', file=sys.stderr):
    """
    Display ay progress bar on stderr like
    
        [===========>           ] 54% of 2300 annotations
    
    :param count: the current number of annotation read
    :type count: int or float 
    :param total: the total number of annotations
    :type totol: int or float
    :param status: the message to display at the right of progress bar
    :type status: str
    :param file: where to write the progress bar (by default sys.stderr). 
    :type file: file like object
    """
    # The MIT License (MIT)
    # Copyright (c) 2016 Vladimir Ignatev
    #
    # Permission is hereby granted, free of charge, to any person obtaining
    # a copy of this software and associated documentation files (the "Software"),
    # to deal in the Software without restriction, including without limitation
    # the rights to use, copy, modify, merge, publish, distribute, sublicense,
    # and/or sell copies of the Software, and to permit persons to whom the Software
    # is furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included
    # in all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
    # FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
    # OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
    # OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    # modified by Bertrand Néron (2017)
    bar_len = 60
    filled_len = int(round(bar_len * count / total))

    percent = count / total
    bar = '=' * filled_len + '>' + ' ' * (bar_len - filled_len - 1)

    file.write('[{bar}] {percent:.1%} ... of {status} annotations\r'.format(bar=bar, percent=percent, status=total))
    file.flush()  # As suggested by Rom Ruben
    # (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)


# pre defined some matloptlib backend in function of output format
non_interactive_backends = {'png': 'AGG',
                            'ps': 'PS',
                            'eps': 'PS',
                            'svg': 'SVG'}


