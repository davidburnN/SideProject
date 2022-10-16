#! /Users/sehuang/opt/anaconda3/bin/python3.8

import readline

class CarryOperation():
    def carryOperationCalculator(self, n1, n2):
        carryOperation = 0
        nl1, nl2 = len(n1), len(n2)
        if nl1<nl2:
            for j in range(nl1):
                if int(n1[nl1-1-j])+int(n2[nl2-1-j])>9:
                    carryOperation += 1
        else:
            for j in range(nl2):
                if int(n1[nl1-1-j])+int(n2[nl2-1-j])>9:
                    carryOperation += 1
        return carryOperation
    
    def OperationCount(self):
        print("Please input couple of positive integers to find out carry operations. At very end please input two 0's in one layer to stop inputing.")
        numbers = []
        carryOperationArray = []
        i = 0

        while True:
            while True:
                try:
                    n1 = input(f'Input the first number in layer {i+1}: ')
                except ValueError:
                    print('Please input a valid number.')
                    continue
                if len(n1)>9:
                    print('Please input an integer less than 10 digits.')
                    continue
                elif int(n1)<0:
                    print('Please input a none negetive integer.')
                    continue
                else:
                    break
            while True:
                try:
                    n2 = input(f'Input the second number in layer {i+1}: ')
                except ValueError:
                    print('Please input a valid number.')
                    continue
                if len(n2)>9:
                    print('Please input an integer less than 10 digits.')
                    continue
                elif int(n2)<0:
                    print('Please input a none negetive integer.')
                    continue
                else:
                    break

            if (int(n1)==0) & (int(n2)==0):
                break
            else:
                numbers.append([n1, n2])
                carryOperationArray.append(CarryOperation().carryOperationCalculator(n1,n2))
                i += 1


        for j in range(i):
            print(numbers[j], f"{carryOperationArray[j]} carry operations")

carryOperation = CarryOperation()
operationCount = carryOperation.OperationCount()
