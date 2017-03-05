#include<iostream>

class Complex{
private:
    double real,
            imag;
public:
    Complex();  //Constructor
    Complex(double);
    Complex(double,double);
    ~Complex(){}; //Destructor  //Same as Complex:~Complex(){}; Below
    
    void print();
    friend Complex operator+(Complex,Complex);
    
    //friend void helloworld();  helloworld function is a friend of the class (Can access private)
}; //number1,number2;

Complex::Complex(){
    real=0.0;
    imag=0.0;
}

Complex::Complex(double re){
    real=re;
    imag=0.0;
}

Complex::Complex(double re, double im){
    real=re;
    imag=im;
}

void Complex::print(){
    std::cout << "(" << real << "," << imag << ")" << std::endl;
}


Complex operator+(Complex a, Complex b){ //Defining '+'
    return Complex(a.real + b.real, a.imag + b.imag);
}
 

int main(int argc, char*argv[]){
    //declare variables
    double real, imag;
    
    //ask for input
    std::cout << "Real: ";
    std::cin >> real;
    std::cout << "Imag: ";
    std::cin >> imag;
    
    Complex number1;
    Complex number2(real);
    Complex number3(real,imag);
    Complex number4 = number2 + number3;
    
    
    /* Will not work since real and imag are private
    std::cout << number3.real << std::endl;
    std::cout << number3.imag << std::endl;
     */
    
    number1.print();
    number2.print();
    number3.print();
    number4.print();
    
    
    
    
    return 0;
}