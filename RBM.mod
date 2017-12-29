/*********************************************
 * OPL 12.7.1.0 Model
 * Author: MattSucks
 * Creation Date: Aug 2, 2017 at 4:20:42 AM
 *********************************************/
 
int qlen = ...;
int rlen = ...;
float weights[1..qlen*rlen] = ...;
 
dvar int+ x[1..qlen*rlen];
 
maximize sum(i in 0..qlen-1, j in 0..rlen-1) weights[i*rlen+j+1]*x[i*rlen+j+1];
 
subject to{
	forall(j in 0..rlen-1){
		sum(i in 0..qlen-1) x[i*rlen+j+1] <= 1;
	}
	forall(i in 0..qlen-1){
		sum(j in 0..rlen-1) x[i*rlen+j+1] <= 1;
	}
}

main {
	thisOplModel.generate();
	cplex.solve();
	var ofile = new IloOplOutputFile("Matching.txt");
	ofile.writeln(thisOplModel.printSolution());
	ofile.close();
}
