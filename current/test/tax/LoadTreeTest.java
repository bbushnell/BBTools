package test.tax;

import parse.Parse;
import tax.TaxTree;

/** Manual taxonomy-loader smoke test; uses configured defaults unless a path is supplied. */
public class LoadTreeTest{
	/** Constructs the stateless manual smoke test. */
	public LoadTreeTest(){}

	/** Loads a tree with optional name indexes and asserts a non-null result.
	 * Run with assertions enabled; this is a loader smoke check, not a content oracle.
	 * @param args Optional tree path, hashNames boolean and hashDot boolean */
	public static void main(String[] args){
		String path="auto";
		boolean hashNames=false, hashDot=false;
		if(args.length>0){path=args[0];}
		if(args.length>1){hashNames=Parse.parseBoolean(args[1]);}
		if(args.length>2){hashDot=Parse.parseBoolean(args[2]);}
//		Timer t=new Timer();
		TaxTree tree=TaxTree.loadTaxTree(path, System.err, hashNames, hashDot);
		assert(tree!=null) : "Tree is null.";
//		t.stopAndPrint();
	}
	
}
