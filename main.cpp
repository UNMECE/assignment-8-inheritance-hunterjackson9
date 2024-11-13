#include <iostream>
#include <cmath>

const double EPSILON_0 = 8.85e-12;  
const double MU_0 = 1.257e-6;       

// Base class Field
class Field {
protected:
    double* value;

public:
    // Default constructor
    Field() {
        value = new double[3];
        value[0] = value[1] = value[2] = 0.0;
    }

    // Constructor with parameters
    Field(double x, double y, double z) {
        value = new double[3];
        value[0] = x;
        value[1] = y;
        value[2] = z;
    }

    // Copy constructor
    Field(const Field& other) {
        value = new double[3];
        for(int i = 0; i < 3; i++) {
            value[i] = other.value[i];
        }
    }

    // Destructor
    virtual ~Field() {
        delete[] value;
    }

    // Print magnitude function
    virtual void printMagnitude() const {
        std::cout << "Components: (" << value[0] << ", " 
                 << value[1] << ", " << value[2] << ")" << std::endl;
    }

    // Assignment operator
    Field& operator=(const Field& other) {
        if(this != &other) {
            for(int i = 0; i < 3; i++) {
                value[i] = other.value[i];
            }
        }
        return *this;
    }
};

// Electric Field class
class Electric_Field : public Field {
private:
    double calculatedField;

public:
    // Default constructor
    Electric_Field() : Field() {
        calculatedField = 0.0;
    }

    // Constructor with parameters
    Electric_Field(double x, double y, double z) : Field(x, y, z) {
        calculatedField = 0.0;
    }

    // Copy constructor
    Electric_Field(const Electric_Field& other) : Field(other) {
        calculatedField = other.calculatedField;
    }

    // Calculate electric field using Gauss's law
    void calculateField(double Q, double r) {
        calculatedField = Q / (4 * M_PI * r * r * EPSILON_0);
    }

    // Get calculated field
    double getCalculatedField() const {
        return calculatedField;
    }

    // Overload + operator
    Electric_Field operator+(const Electric_Field& other) const {
        return Electric_Field(
            value[0] + other.value[0],
            value[1] + other.value[1],
            value[2] + other.value[2]
        );
    }

    // Friend function for << operator overload
    friend std::ostream& operator<<(std::ostream& os, const Electric_Field& ef);
};

// Magnetic Field class
class Magnetic_Field : public Field {
private:
    double calculatedField;

public:
    // Default constructor
    Magnetic_Field() : Field() {
        calculatedField = 0.0;
    }

    // Constructor with parameters
    Magnetic_Field(double x, double y, double z) : Field(x, y, z) {
        calculatedField = 0.0;
    }

    // Copy constructor
    Magnetic_Field(const Magnetic_Field& other) : Field(other) {
        calculatedField = other.calculatedField;
    }

    // Calculate magnetic field using Ampere's law
    void calculateField(double I, double r) {
        calculatedField = (MU_0 * I) / (2 * M_PI * r);
    }

    // Get calculated field
    double getCalculatedField() const {
        return calculatedField;
    }

    // Overload + operator
    Magnetic_Field operator+(const Magnetic_Field& other) const {
        return Magnetic_Field(
            value[0] + other.value[0],
            value[1] + other.value[1],
            value[2] + other.value[2]
        );
    }

    // Friend function for << operator overload
    friend std::ostream& operator<<(std::ostream& os, const Magnetic_Field& mf);
};

// Overload << operator for Electric_Field
std::ostream& operator<<(std::ostream& os, const Electric_Field& ef) {
    os << "Electric Field Components: (" << ef.value[0] << ", " 
       << ef.value[1] << ", " << ef.value[2] << ")";
    return os;
}

// Overload << operator for Magnetic_Field
std::ostream& operator<<(std::ostream& os, const Magnetic_Field& mf) {
    os << "Magnetic Field Components: (" << mf.value[0] << ", " 
       << mf.value[1] << ", " << mf.value[2] << ")";
    return os;
}

int main() {
    // Create electric field objects
    Electric_Field e1(1.0, 2.0, 3.0);
    Electric_Field e2(4.0, 5.0, 6.0);
    
    // Create magnetic field objects
    Magnetic_Field m1(0.1, 0.2, 0.3);
    Magnetic_Field m2(0.4, 0.5, 0.6);

    // Print magnitudes using base class function
    std::cout << "Electric Field 1:" << std::endl;
    e1.printMagnitude();
    std::cout << "Electric Field 2:" << std::endl;
    e2.printMagnitude();
    
    std::cout << "\nMagnetic Field 1:" << std::endl;
    m1.printMagnitude();
    std::cout << "Magnetic Field 2:" << std::endl;
    m2.printMagnitude();

    // Calculate fields using Gauss's and Ampere's laws
    double charge = 1.0e-6;    // 1 μC
    double current = 2.0;      // 2 A
    double distance = 0.1;     // 0.1 m

    e1.calculateField(charge, distance);
    m1.calculateField(current, distance);

    std::cout << "\nCalculated Electric Field at r = " << distance << "m: " 
              << e1.getCalculatedField() << " N/C" << std::endl;
    std::cout << "Calculated Magnetic Field at r = " << distance << "m: " 
              << m1.getCalculatedField() << " T" << std::endl;

    // Demonstrate operator overloading
    Electric_Field e3 = e1 + e2;
    Magnetic_Field m3 = m1 + m2;

    std::cout << "\nSum of Electric Fields:" << std::endl;
    std::cout << e3 << std::endl;

    std::cout << "\nSum of Magnetic Fields:" << std::endl;
    std::cout << m3 << std::endl;

    return 0;
}
