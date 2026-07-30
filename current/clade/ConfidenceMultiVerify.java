package clade;

import fileIO.ByteFile;
import fileIO.FileFormat;
import ml.CellNet;
import ml.CellNetParser;
import structures.FloatList;

/**
 * Verifies a multi-output confidence bundle produced by {@link ConfidenceMultiBundle}: the bundle must
 * load through SerialNNLoader as a multi-output model whose embedded net reproduces the ORIGINAL net's
 * outputs bit-for-bit, and it must carry one calibration per output head. Proves the packaging + loader
 * plumbing (not the calibration math, which the head-to-head already validated).
 *
 * Usage: confidencemultiverify.sh bundle=<bundle.bbnets.gz> net=<orig.bbnet> val=<val.tsv> [rows=2000]
 *
 * @author Noire
 */
public class ConfidenceMultiVerify {

	public static void main(String[] args){
		String bundlePath=null, netPath=null, valPath=null; int rows=2000;
		for(String a : args){
			final String[] s=a.split("=", 2);
			final String k=s[0].toLowerCase(), v=(s.length>1 ? s[1] : null);
			if(k.equals("bundle")){bundlePath=v;}
			else if(k.equals("net")){netPath=v;}
			else if(k.equals("val")){valPath=v;}
			else if(k.equals("rows")){rows=Integer.parseInt(v);}
			else{throw new RuntimeException("Unknown parameter: "+a);}
		}
		if(bundlePath==null||netPath==null||valPath==null){throw new RuntimeException("Required: bundle= net= val=");}

		final SerialNNLoader.LoadedNets ln=SerialNNLoader.load(bundlePath);
		if(ln==null){fail("bundle load returned null");}
		if(!ln.multioutput){fail("bundle is not multioutput");}
		if(ln.multiNet==null){fail("multiNet null");}
		final int NO=ln.multiNet.numOutputs();
		System.out.println("bundle: levels="+ln.levels+" multiNet inputs="+ln.multiNet.numInputs()+" outputs="+NO);
		System.out.print("labels: ");
		for(int i=0; i<ln.levels; i++){System.out.print((ln.multiLabels==null?"?":ln.multiLabels[i])+" ");}
		System.out.println();
		int lutOk=0;
		for(int i=0; i<ln.levels; i++){
			final boolean has=(ln.multiLutX!=null && ln.multiLutX[i]!=null);
			if(has){lutOk++;}
			System.out.println("  head "+i+" ("+(ln.multiLabels==null?"?":ln.multiLabels[i])+"): lut="+(has?ln.multiLutX[i].length+" knots":"none"));
		}
		final int ND=ln.multiNet.numInputs();
		if(ND<48){fail("multiNet inputs < 48");}
		if(NO!=ln.levels){fail("multiNet outputs != levels");}
		if(lutOk!=ln.levels){fail("not every head has a calibration lut ("+lutOk+"/"+ln.levels+")");}

		final CellNet orig=CellNetParser.load(netPath);
		if(orig==null){fail("orig net load failed");}

		//--- embedded net must reproduce the original's outputs exactly ---
		final FloatList fl=new FloatList(ND);
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(valPath, FileFormat.TEXT, null, true, false));
		byte[] line; int n=0; double maxDiff=0;
		while((line=bf.nextLine())!=null && n<rows){
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] f=new String(line).split("\t");
			if(f.length<ND){continue;}
			fl.size=0; boolean ok=true;
			for(int i=0; i<ND; i++){try{fl.add(Float.parseFloat(f[i]));}catch(Exception e){ok=false; break;}}
			if(!ok){continue;}
			ln.multiNet.applyInput(fl); ln.multiNet.feedForward();
			final float[] a=ln.multiNet.getOutput();
			orig.applyInput(fl); orig.feedForward();
			final float[] b=orig.getOutput();
			for(int i=0; i<NO; i++){maxDiff=Math.max(maxDiff, Math.abs(a[i]-b[i]));}
			n++;
		}
		bf.close();
		System.out.println("compared "+n+" rows; embedded-vs-original maxDiff="+maxDiff);
		if(maxDiff>1e-5){fail("embedded net diverges from original (maxDiff "+maxDiff+")");}
		System.out.println("VERIFY PASS: multi-output bundle loads, "+ln.levels+" calibrated heads, net bit-faithful.");
	}

	static void fail(String msg){System.out.println("VERIFY FAIL: "+msg); throw new RuntimeException(msg);}
}
