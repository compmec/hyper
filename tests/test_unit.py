
import sys
import traceback


class TU:
    """
    Test unity class
    """
    SUCESS = True
    FAILURE = False
    FOLDER_SRC = "../src/"
    FOLDER_GEO = "../geo/"
    FOLDER_MSH = "../msh/"
    VERBOSE = True
    NUM = 0
    DEN = 0
    FUNCTIONS = []

    def __init__(self):
        self.result = TU.SUCESS
        self.names = []
        self.numbers = []
        self.functions = []
        if sys.platform == "linux" or sys.platform == "linux2":
            # linux
            pass
        elif sys.platform == "darwin":
            # OS X
            pass
        elif sys.platform == "win32":
            # Windows...
            # print("windows!")
            pass

    @staticmethod
    def print_result(msg, result):
        print(msg + ": " + str("SUCESS" if result == TU.SUCESS else "FAILURE"))

    @staticmethod
    def print_final_result(msg, result):
        TU.print_result(msg, result)

    @staticmethod
    def include_path(name):
        if name == "src":
            sys.path.append(TU.FOLDER_SRC)
        elif name == "geo":
            sys.path.append(TU.FOLDER_GEO)
        elif name == "msh":
            sys.path.append(TU.FOLDER_MSH)

    @staticmethod
    def get_numberoftests(function):
        number = 0
        while 1:
            try:
                function(number)
            except TError:
                return number
            except:
                pass
            finally:
                number += 1
            if number > 100:
                raise Exception("Could not get the exemples")
        return number

    def get_classname(self):
        return str(type(self)).split("'")[1].split(".")[-1]

    def __test_one_function(self, i):
        test_function = self.functions[i]
        try:
            test_function()
        except Exception as e:  # Doesn't matter the exception, it's an error
            if self.VERBOSE:
                traceback.print_exc()
                # print("")
                # print("")
                # print("#" * 50)
                # print("# Error given by:")
                # print("#     " + str(e))
                # print("#" * 50)
            self.result = self.FAILURE
        if self.VERBOSE:
            self.NUM += 1
            index = "[" + str(self.NUM) + "/" + str(self.DEN) + "]"
            name = test_function.__name__
            if name.startswith("__"):
                name = name[2:]
            self.print_result("    " + index + " " + name, self.result)

    def run(self):

        classname = self.get_classname()
        setUp_name = "_" + str(classname) + "__setUp"
        tearDown_name = "_" + str(classname) + "__tearDown"

        if setUp_name in dir(self):
            self.DEN += 1
            setUp_func = getattr(self, setUp_name)
            self.functions.append(setUp_func)

        for function_name in self.FUNCTIONS:
            self.names.append(function_name)
            func = getattr(self, function_name)
            self.functions.append(func)

        self.DEN += len(self.FUNCTIONS)
        if tearDown_name in dir(self):
            self.DEN += 1
            tearDown_func = getattr(self, tearDown_name)
            self.functions.append(tearDown_func)

        if self.VERBOSE:
            print(classname + ": Routine testing")
        for i in range(self.DEN):
            self.__test_one_function(i)
            if self.result == self.FAILURE:
                break
        self.print_final_result(classname, self.result)
        return self.result


class TError(BaseException):
    """
    Exception for verify quantity of tests
    """
    # pass

    def __init__(self, message):
        self.message = message
        super(TError, self).__init__()

    def __str__(self):
        return self.message
