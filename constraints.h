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
		: constraintType(_constraintType),
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

		MatrixXd constGradient = MatrixXd::Zero(4, 3);
		RowVector3d r1 = currVertexPositions.row(0) - currCOMPositions.row(0);
		RowVector3d r2 = currVertexPositions.row(1) - currCOMPositions.row(1);
		RowVector3d n;

		switch (constraintType)
		{
			case ConstraintType::COLLISION:
			{
				n = refVector;
			}
			break;
			case ConstraintType::DISTANCE:
			{
				RowVector3d v1 = currVertexPositions.row(0);
				RowVector3d v2 = currVertexPositions.row(1);
				auto x = v2 - v1;
				double dist = sqrt(x.dot(x));
				n = x / dist;
			}
			break;
		}
		constGradient.row(0) = -n;
		constGradient.row(1) = -r1.cross(n);
		constGradient.row(2) = n;
		constGradient.row(3) = r2.cross(n);

		constGradient = constGradient.transpose().eval();
		constGradient.resize(1, 12);

		/**************
		 TODO: write velocity correction procedure:
		 1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true

		 2. Otherwise, correct linear and angular velocities as learnt in class.

		 Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
		 ***************/

		 // 1x12 vector containing v1, omega 1, v2, omega 2 
		MatrixXd v(4, 3);
		// v1
		v.row(0) = currCOMVelocities.row(0);
		// omega 1
		v.row(1) = currAngularVelocities.row(0);
		// v2
		v.row(2) = currCOMVelocities.row(1);
		// omega 2
		v.row(3) = currAngularVelocities.row(1);

		v = v.transpose().eval();
		v.resize(12, 1);

		double Jv = (constGradient * v)(0, 0);

		if ((constraintEqualityType == ConstraintEqualityType::EQUALITY && std::abs(Jv) <= tolerance)
			|| (constraintEqualityType == ConstraintEqualityType::INEQUALITY && Jv >= 0))
		{
			correctedCOMVelocities = currCOMVelocities;
			correctedAngularVelocities = currAngularVelocities;
			return true;
		}

		// TEMP
		RowVector3d r1n = r1.cross(n);
		RowVector3d r2n = r2.cross(n);
		double inertia1 = r1n * invInertiaTensor1 * r1n.transpose();
		double inertia2 = r2n * invInertiaTensor2 * r2n.transpose();
		RowVector3d ang1 = currAngularVelocities.row(0);
		RowVector3d ang2 = currAngularVelocities.row(1);

		RowVector3d v1Min = currCOMVelocities.row(0) + ang1.cross(r1);
		RowVector3d v2Min = currCOMVelocities.row(1) + ang2.cross(r2);

		double nomP1 =  -(1 + CRCoeff) * (v2Min - v1Min).dot(n);
		double denomP1 = (invMass1 + invMass2 + inertia1 + inertia2);

		double nomP2 = -(1 + CRCoeff) * Jv;
		double denomP2 = (constGradient * invMassMatrix * constGradient.transpose())(0, 0);
		double lambda = nomP2 / denomP2;

		// [com1 ang1 com2 ang2]
		// resize to 4x3
		MatrixXd corrections;
		corrections = (invMassMatrix * constGradient.transpose()) * lambda;
		corrections.resize(3, 4);
		corrections = corrections.transpose().eval();

		MatrixXd deltaCOM(2, 3);
		MatrixXd deltaANG(2, 3);
		deltaCOM.row(0) = corrections.row(0);
		deltaCOM.row(1) = corrections.row(2);
		deltaANG.row(0) = corrections.row(1);
		deltaANG.row(1) = corrections.row(3);

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
		auto x = v2 - v1;
		double dist = sqrt(x.dot(x));

		double Cp;

		MatrixXd J(2, 3);
		RowVector3d n;
		switch (constraintType)
		{
			case ConstraintType::COLLISION:
			{
				n = refVector;
				Cp = -dist;
			}
			break;
			case ConstraintType::DISTANCE:
			{
				// The current distance (dist), and the desired distance (refValue), must be equal to zero.
				n = x / dist;
				Cp = x.dot(n) - refValue;
			}
			break;
		}
		J.row(0) = -n;
		J.row(1) = n;
		J = J.transpose().eval();
		J.resize(1, 6);

		if ((constraintEqualityType == ConstraintEqualityType::EQUALITY && std::abs(Cp) <= tolerance)
			|| (constraintEqualityType == ConstraintEqualityType::INEQUALITY && Cp >= 0))
		{
			correctedCOMPositions = currCOMPositions;
			return true;
		}

		double frac = (J * invMassMatrix * J.transpose())(0, 0);
		double lambda = Cp / frac;

		// 1x6
		correctedCOMPositions = -invMassMatrix * J.transpose() * lambda;
		// Resize to 2x3 
		correctedCOMPositions.resize(3, 2);
		correctedCOMPositions = correctedCOMPositions.transpose().eval();

		//RowVector3d temp = correctedCOMPositions.row(0);
		correctedCOMPositions += currCOMPositions;

		return false;
	}
};
#endif /* constraints_h */
