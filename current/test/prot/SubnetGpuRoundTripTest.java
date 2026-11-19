package test.prot;

import java.util.Locale;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;

/**
 * Compares exported Java predictions with saved GPU predictions or a second
 * native network. Reference-network mode requires exact output bits and also
 * checks that compression/metadata changes preserve all declared output heads.
 * Raw validation inputs are read directly, including their declared target columns;
 * Only input columns enter CellNet. Declared training targets may differ from the
 * network's output count because live error heads need no extra stored labels.
 * No normalization is recomputed here.
 * @author Yoimiya
 */
public final class SubnetGpuRoundTripTest {

	/** Requires dimensions and either saved predictions or an exact native reference. */
	public static void main(String[] args){
		String netFile=null,vectors=null,predictions=null,referenceFile=null;
		int rows=1000;
		double tolerance=0.00002;
		for(String arg : new PreParser(args,null,false).args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("net")){netFile=value;}
			else if(key.equals("vectors")){vectors=value;}
			else if(key.equals("predictions")){predictions=value;}
			else if(key.equals("reference")){referenceFile=value;}
			else if(key.equals("rows")){rows=Integer.parseInt(value);}
			else if(key.equals("tolerance")){tolerance=Double.parseDouble(value);}
			// Exported TSV values may use exponents; PreParser does not consume this common flag.
			else if(key.equals("forcejavaparsedouble") && Parser.parseCommonStatic(arg,key,value)){}
			else{throw new IllegalArgumentException("Unknown argument: "+key);}
		}
		if(netFile==null || vectors==null || (predictions==null)==(referenceFile==null) || rows<1 ||
				!Double.isFinite(tolerance) || tolerance<0 || (predictions!=null && tolerance==0)){
			throw new IllegalArgumentException("Required: net= vectors= and exactly one of predictions= or reference=; positive rows and valid tolerance");
		}
		final CellNet net=CellNetParser.loadInferenceFromLines(ByteFile.toLines(netFile));
		if(net==null){throw new IllegalArgumentException("Could not load network: "+netFile);}
		assert(net.numInputs()>0 && net.numOutputs()>0) : "CellNet inference requires nonempty input and output layers";
		final CellNet reference=referenceFile==null ? null : CellNetParser.loadInferenceFromLines(ByteFile.toLines(referenceFile));
		if(referenceFile!=null && reference==null){throw new IllegalArgumentException("Could not load reference network: "+referenceFile);}
		final double maximum=compare(net,vectors,predictions,rows,tolerance,reference);
		System.out.println("status\t"+(referenceFile==null ? "JAVA_GPU_PREDICTION_PASS" : "JAVA_JAVA_EXACT_PREDICTION_PASS"));
		System.out.println("rows\t"+rows);
		System.out.println("inputs\t"+net.numInputs());
		System.out.println("outputs\t"+net.numOutputs());
		System.out.println("max_absolute_difference\t"+maximum);
		System.out.println("tolerance\t"+tolerance);
	}

	/** Streams aligned raw/GPU rows and fails on truncation, extra predictions, or mismatch. */
	private static double compare(CellNet net,String vectors,String predictions,int rows,double tolerance){
		return compare(net,vectors,predictions,rows,tolerance,null);
	}

	/** Checks exact native equality with an optional additional saved-prediction bound. */
	public static double compare(CellNet net,String vectors,String predictions,int rows,double tolerance,CellNet loose){
		if(predictions==null && loose==null){throw new IllegalArgumentException("A prediction file or reference network is required");}
		if(loose!=null && (loose.numInputs()!=net.numInputs() || loose.numOutputs()!=net.numOutputs())){
			throw new IllegalArgumentException("Loose and bundled network dimensions differ");
		}
		final ByteFile vf=ByteFile.makeByteFile(vectors,true);
		ByteFile pf=null;
		try{
			pf=predictions==null ? null : ByteFile.makeByteFile(predictions,true);
			final LineParser1 v=new LineParser1((byte)'\t'),p=new LineParser1((byte)'\t');
			final byte[] header=vf.nextLine();
			if(header==null){throw new IllegalArgumentException("Empty vector file");}
			v.set(header);
			if(v.terms()!=4 || !v.termEquals("#dims",0) || v.parseInt(1)!=net.numInputs()
					|| v.parseInt(2)<0 || v.parseInt(3)!=0){
				throw new IllegalArgumentException("Vector dimensions do not match exported network inputs");
			}
			final int targetColumns=v.parseInt(2);
			// Live residual and classification heads derive targets from the two stored value labels.
			if(targetColumns!=net.numOutputs() && !(targetColumns==2 &&
					(net.numOutputs()==4 || net.numOutputs()==6))){
				throw new IllegalArgumentException("Stored targets must match network outputs, or use two value labels for four or six live outputs");
			}
			final float[] input=new float[net.numInputs()];
			double maximum=0;
			for(int row=0; row<rows; row++){
				final byte[] vl=nextData(vf),pl=pf==null ? null : nextData(pf);
				if(vl==null || (pf!=null && pl==null)){throw new IllegalArgumentException("Missing raw/reference row "+row);}
				v.set(vl); if(pf!=null){p.set(pl);}
				if(v.terms()!=net.numInputs()+targetColumns || (pf!=null && (p.terms()!=net.numOutputs()+1
						|| p.parseInt(0)!=row))){
					throw new IllegalArgumentException("Input/prediction shape or row identity mismatch at "+row);
				}
				for(int i=0; i<input.length; i++){
					input[i]=v.parseFloat(i);
					if(!Float.isFinite(input[i])){throw new IllegalArgumentException("Nonfinite input at row "+row+" column "+i);}
				}
				net.applyInput(input);
				net.feedForward();
				if(loose!=null){loose.applyInput(input); loose.feedForward();}
				for(int i=0; i<net.numOutputs(); i++){
					if(loose!=null && Float.floatToIntBits(net.getOutput(i))!=Float.floatToIntBits(loose.getOutput(i))){
						throw new IllegalStateException("Loose/bundled output bits differ at row="+row+" output="+i);
					}
					final double expected=pf==null ? loose.getOutput(i) : p.parseDouble(i+1),actual=net.getOutput(i);
					final double error=Math.abs(actual-expected);
					if(!Double.isFinite(error) || error>tolerance){
						throw new IllegalStateException("GPU/Java prediction mismatch row="+row+" output="+i+
							" expected="+expected+" actual="+actual+" absolute_error="+error+" tolerance="+tolerance);
					}
					maximum=Math.max(maximum,error);
				}
			}
			if(pf!=null && nextData(pf)!=null){throw new IllegalArgumentException("GPU predictions exceed requested row count "+rows);}
			assert(maximum<=tolerance) : "Every compared output must satisfy the absolute-error bound";
			return maximum;
		}finally{
			final boolean vError=vf.close(),pError=pf!=null && pf.close();
			if(vError || pError){throw new IllegalStateException("I/O error closing prediction comparison inputs");}
		}
	}

	/** Skips only empty/comment lines; data order must remain unchanged. */
	private static byte[] nextData(ByteFile bf){
		assert(bf!=null) : "Prediction comparison must open its input before reading";
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]!='#'){return line;}
		}
		return null;
	}
}
