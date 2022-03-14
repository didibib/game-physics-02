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
		RowVectorXd constGradient(12);

		/**************
		 TODO: write velocity correction procedure:
		 1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true

		 2. Otherwise, correct linear and angular velocities as learnt in class.

		 Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
		 ***************/
		// 1x12 vector containing v1, omega 1, v2, omega 2 
		RowVectorXd v;
		// v1
		v << currCOMVelocities.row(0);
		// omega 1
		v << currAngularVelocities.row(0);
		// v2
		v << currCOMVelocities.row(1);
		// omega 2
		v << currAngularVelocities.row(1);	

		//Stub code: remove upon implementation
		correctedCOMVelocities = currCOMVelocities;
		correctedAngularVelocities = currAngularVelocities;
		return true;
		//end of stub code
	}

	//projects the position unto the constraint
	//returns true if constraint was already valid with "currPositions"
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
		RowVectorXd J;
		J << refVector;
		J << -refVector;

		double lambda = Cp / (J * invMassMatrix * J.transpose());

		// 1x6
		correctedCOMPositions = invMassMatrix * J.transpose() * lambda;
		// Resize to 2x3 
		correctedCOMPositions.resize(2, 3);

		return true;
	}
};



#endif /* constraints_h */
