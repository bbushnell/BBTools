package cardinality;

/**
 * FarmEval: scores a train.sh-trained CellNet on the farm's held-out family
 * vector files (written by ComplexityFarm with a dump path).
 * Target denormalization: y = yNorm*13 - 13 (log2 complexity), cHat = 2^y.
 * Reports mean |cHat - c| / c per family — directly comparable to the
 * quad%/bigNN% columns in farm output.
 *
 * Usage: java cardinality.FarmEval <net.bbnet> <evalFile1.tsv> [evalFile2.tsv ...]
 *
 * @author Amber
 * @date July 2026
 */
public class FarmEval {

	public static void main(String[] args) throws Exception{
		final ml.CellNet net=ml.CellNetParser.load(args[0]);
		if(net==null){throw new RuntimeException("cannot load net "+args[0]);}
		System.out.println("#family\tsamples\tcellnet_cErr%");
		double tot=0;
		int fams=0;
		for(int a=1; a<args.length; a++){
			double errSum=0;
			int n=0;
			for(String line : java.nio.file.Files.readAllLines(java.nio.file.Paths.get(args[a]))){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] f=line.split("\t");
				final float[] x=new float[f.length-1];
				for(int i=0; i<x.length; i++){x[i]=Float.parseFloat(f[i]);}
				final double yTrue=Double.parseDouble(f[f.length-1])*13-13;
				net.applyInput(x);
				final double yHat=net.feedForward()*13-13;
				final double c=Math.pow(2.0, yTrue), cHat=Math.pow(2.0, yHat);
				errSum+=Math.abs(cHat-c)/c;
				n++;
			}
			final String name=args[a].replaceAll(".*eval_", "").replaceAll("\\.tsv$", "");
			final double e=100*errSum/Math.max(1, n);
			System.out.println(String.format("%s\t%d\t%.2f", name, n, e));
			tot+=e; fams++;
		}
		System.out.println(String.format("#MEAN\t-\t%.2f", tot/Math.max(1, fams)));
	}
}
