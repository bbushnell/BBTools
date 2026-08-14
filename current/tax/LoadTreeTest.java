package tax;

import parse.Parse;

public class LoadTreeTest{

	public static void main(String[] args) {
		String path="auto";
		boolean hashNames=false, hashDot=false;
		if(args.length>0) {path=args[0];}
		if(args.length>1) {hashNames=Parse.parseBoolean(args[1]);}
		if(args.length>2) {hashDot=Parse.parseBoolean(args[2]);}
//		Timer t=new Timer();
		TaxTree tree=TaxTree.loadTaxTree(path, System.err, hashNames, hashDot);
		assert(tree!=null) : "Tree is null.";
//		t.stopAndPrint();
	}
	
}
