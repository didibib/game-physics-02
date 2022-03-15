#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum class ConstraintType { DISTANCE, COLLISION } ConstraintType;   //You can expand it for more constraints
typedef enum class ConstraintEqualityType { EQUALITY, INEQUALITY } ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint
{
public:

	int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
	int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
	double invMass1, invMass2;       //inverse masses of two bodies
	double refValue;                //Reference values to use in the constraint, when needed (like distance)
	RowVector3d refVector;          //Reference vector when needed (like vector)
	double CRCoeff;                 //extra velocity bias
	ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
	ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality

	Constraint(const ConstraintType _constraintType, 
			   const ConstraintEqualityType _constraintEqualityType,
			   const int& _m1,
			   const int& _v1,
			   const int& _m2,
			   const int& _v2,
			   const double& _invMass1,
			   const double& _invMass2,
			   const RowVector3d& _refVector,
			   const double& _refValue,
			   const double& _CRCoeff)
		:	constraintType(_constraintType),
			constraintEqualityType(_constraintEqualityType),
			m1(_m1), v1(_v1), m2(_m2), v2(_v2),
			invMass1(_invMass1),
			invMass2(_invMass2),
			refValue(_refValue),
			CRCoeff(_CRCoeff)
	{
		refVector = _refVector;
	}

	~Constraint() {}

	//computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
	//The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
	//returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
	//currCOMPositions is a 2x3 matrix, where each row is per one of the sides of the constraints; the rest of the relevant variables are similar, and so should the outputs be resized.
	bool resolveVelocityConstraint(const MatrixXd& currCOMPositions,
								   const MatrixXd& currVertexPositions,
								   const MatrixXd& currCOMVelocities,
								   const MatrixXd& currAngularVelocities,
								   const Matrix3d& invInertiaTensor1,
								   const Matrix3d& invInertiaTensor2,
										 MatrixXd& correctedCOMVelocities,
									     MatrixXd& correctedAngularVelocities,
										 double    tolerance)
	{
		MatrixXd invMassMatrix = MatrixXd::Zero(12, 12);
		invMassMatrix(0, 0) = invMass1;
		invMassMatrix(1, 1) = invMass1;
		invMassMatrix(2, 2) = invMass1;

		invMassMatrix(6, 6) = invMass2;
		invMassMatrix(7, 7) = invMass2;
		invMassMatrix(8, 8) = invMass2;

		invMassMatrix.block(3, 3, 3, 3) = invInertiaTensor1;
		invMassMatrix.block(9, 9, 3, 3) = invInertiaTensor2;

		MatrixXd constGradient(4,3);

		RowVector3d r1 = currVertexPositions.row(0) - currCOMPositions.row(0);
		RowVector3d r2 = currVertexPositions.row(1) - currCOMPositions.row(1);

		// normal
		constGradient.row(0) = refVector;
		constGradient.row(1) = refVector.cross(r1);
		constGradient.row(2) = -refVector;
		constGradient.row(3) = -refVector.cross(r2);

		constGradient.resize(1, 12);

		/**************
		 TODO: write velocity correction procedure:
		 1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true

		 2. Otherwise, correct linear and angular velocities as learnt in class.

		 Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
		 ***************/

		// 1x12 vector containing v1, omega 1, v2, omega 2 
		MatrixXd v (4,3);
		// v1
		v.row(0) = currCOMVelocities.row(0);
		// omega 1
		v.row(1) = currAngularVelocities.row(0);
		// v2
		v.row(2) = currCOMVelocities.row(1);
		// omega 2
		v.row(3) = currAngularVelocities.row(1);	
		v.resize(12, 1);
		double Jv = (constGradient * v)(0,0);

		if ((constraintEqualityType == ConstraintEqualityType::EQUALITY && std::abs(Jv) <= tolerance)
			|| (constraintEqualityType == ConstraintEqualityType::INEQUALITY && Jv >= 0))
		{
			correctedCOMVelocities = currCOMVelocities;
			correctedAngularVelocities = currAngularVelocities;
			return true;
		}

		double frac = (constGradient * invMassMatrix * constGradient.transpose())(0,0);
		double lambda = -(1 + CRCoeff) * Jv / frac;

		// [com1 ang1 com2 ang2]
		RowVector3d corrections = (invMassMatrix * constGradient.transpose()) * lambda;
		corrections.resize(4, 3);

		MatrixXd deltaCOM(2, 3);
		MatrixXd deltaANG(2, 3);
		deltaCOM << corrections.row(0);
		deltaCOM << corrections.row(2);
		deltaANG << corrections.row(1);
		deltaANG << corrections.row(3);

		correctedCOMVelocities = currCOMVelocities + deltaCOM;
		correctedAngularVelocities = currAngularVelocities + deltaANG;
		return false;
	}

	// Projects the position unto the constraint
	// Returns true if constraint was already valid with "currPositions"
	bool resolvePositionConstraint(const MatrixXd& currCOMPositions, 
								   const MatrixXd& currConstPositions, 
								         MatrixXd& correctedCOMPositions, 
										 double    tolerance)
	{
		MatrixXd invMassMatrix = MatrixXd::Zero(6, 6);
		invMassMatrix.diagonal() << invMass1, invMass1, invMass1, invMass2, invMass2, invMass2;
		
		/**************
		 TODO: write position correction procedure:
		 1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones and return true

		 2. Otherwise, correct COM position as learnt in class. Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.

		 Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
		***************/

		// Points on the mesh, on which the constraint operates 
		RowVector3d v1 = currConstPositions.row(0);
		RowVector3d v2 = currConstPositions.row(1);
		auto tmp = v2 - v1;
		double dist = tmp.dot(tmp);

		double Cp;
		switch(constraintType)
		{
			case ConstraintType::COLLISION:
			{
				Cp = -dist;
			}			  
			break;
			case ConstraintType::DISTANCE:
			{
				// The current distance (dist), and the desired distance (refValue), must be equal to zero.
				Cp = dist - refValue;
			}
			break;
		}

		if ((constraintEqualityType == ConstraintEqualityType::EQUALITY && std::abs(Cp) <= tolerance)
			|| (constraintEqualityType == ConstraintEqualityType::INEQUALITY && Cp >= 0))
		{
			correctedCOMPositions = currCOMPositions;
			return true;
		}

		// 1x6
		MatrixXd J(2,3);
		J.row(0) = refVector;
		J.row(1) = -refVector;
		J.resize(1, 6);
		
		double frac = (J * invMassMatrix * J.transpose())(0,0);
		double lambda = Cp / frac;

		// 1x6
		correctedCOMPositions = invMassMatrix * J.transpose() * lambda;
		// Resize to 2x3 
		correctedCOMPositions.resize(2, 3);
		correctedCOMPositions += currCOMPositions;

		return false;
	}
};



#endif /* constraints_h */
