package clade;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import ml.CellNet;
import ml.CellNetParser;
import structures.ByteBuilder;
import structures.FloatList;

/**
 * Packages a single multi-output confidence net (RegressionTrainer, 48 inputs -> 8 per-level outputs)
 * into a .bbnets bundle that {@link CladeConfidence} loads as a ONE-net replacement for the 8/9-net
 * per-level bundle. One forward pass, eight heads.
 *
 * The bundle embeds the net verbatim (its own .bbnet lines, re-parsed by CellNetParser) and adds a
 * per-output ISOTONIC calibration table (PAV, output_i -> P(correct at level i)) fit on the training
 * data -- the same calibration the per-level bundle carries, just eight tables against one net.
 *
 * Bundle format:
 *   #levels 8
 *   #multioutput 1
 *   ##network
 *   &lt;embedded .bbnet lines&gt;
 *   #label species
 *   #lut &lt;n&gt;
 *   &lt;x&gt;\t&lt;y&gt;   (n rows)
 *   ... (8 #label/#lut pairs, in output order species..domain)
 *
 * Usage: confidencemultibundle.sh net=<8out.bbnet> fit=<reg_train_ml8.tsv> out=<bundle.bbnets.gz>
 *        [labels=species,genus,family,order,class,phylum,kingdom,domain] [fitmax=300000] [knots=64]
 *
 * @author Noire
 */
public class ConfidenceMultiBundle {

	static final String[] DEF_LABELS={"species","genus","family","order","class","phylum","kingdom","domain"};

	public static void main(String[] args){
		String netPath=null, fitPath=null, outPath=null;
		String[] labels=DEF_LABELS;
		int fitmax=300000, knots=64;
		for(String a : args){
			final String[] s=a.split("=", 2);
			final String k=s[0].toLowerCase(), v=(s.length>1 ? s[1] : null);
			if(k.equals("net")){netPath=v;}
			else if(k.equals("fit") || k.equals("in")){fitPath=v;}
			else if(k.equals("out")){outPath=v;}
			else if(k.equals("labels")){labels=v.split(",");}
			else if(k.equals("fitmax")){fitmax=Integer.parseInt(v);}
			else if(k.equals("knots")){knots=Integer.parseInt(v);}
			else{throw new RuntimeException("Unknown parameter: "+a);}
		}
		if(netPath==null||fitPath==null||outPath==null){throw new RuntimeException("Required: net= fit= out=");}
		final int NL=labels.length;

		//--- load the net twice: as raw lines (to embed) and as a CellNet (to score for calibration) ---
		final ArrayList<byte[]> netLines=ByteFile.toLines(netPath);
		if(netLines==null || netLines.isEmpty()){throw new RuntimeException("empty net file: "+netPath);}
		final CellNet net=CellNetParser.loadFromLines(netLines);
		if(net==null){throw new RuntimeException("net parse failed: "+netPath);}
		if(net.numOutputs()!=NL){
			throw new RuntimeException("net has "+net.numOutputs()+" outputs but "+NL+" labels given");
		}
		final int ND=net.numInputs();   //48, or 49 with the ranking dim

		//--- collect per-output (score, label) over the fit sample; fit rows are ND feats + NL booleans ---
		final FloatList fl=new FloatList(ND);
		final ArrayList<float[]> scoreRows=new ArrayList<float[]>();   //float[NL] net outputs per fit row
		final ArrayList<byte[]> labRows=new ArrayList<byte[]>();       //byte[NL] labels per fit row
		{
			final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(fitPath, FileFormat.TEXT, null, true, false));
			byte[] line;
			while((line=bf.nextLine())!=null && scoreRows.size()<fitmax){
				if(line.length==0 || line[0]=='#'){continue;}
				final String[] f=new String(line).split("\t");
				if(f.length!=ND+NL){continue;}
				fl.size=0;
				boolean ok=true;
				for(int i=0; i<ND; i++){try{fl.add(Float.parseFloat(f[i]));}catch(Exception e){ok=false; break;}}
				if(!ok){continue;}
				net.applyInput(fl); net.feedForward();
				final float[] sc=net.getOutput();
				final byte[] lab=new byte[NL];
				try{for(int L=0; L<NL; L++){lab[L]=(byte)(Float.parseFloat(f[ND+L])>=0.5f?1:0);}}
				catch(Exception e){continue;}
				scoreRows.add(sc); labRows.add(lab);
			}
			bf.close();
		}
		final int nf=scoreRows.size();
		System.err.println("fit rows="+nf+"  outputs="+NL);

		//--- write the bundle ---
		final ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(outPath, FileFormat.TEXT, null, true, true, false, false));
		bsw.start();
		bsw.println(new ByteBuilder("#levels "+NL));         //space-delimited, matching SerialNNLoader
		bsw.println(new ByteBuilder("#multioutput 1"));
		bsw.println(new ByteBuilder("##network 0 -1"));       //single net; level/bin placeholders
		for(byte[] nl : netLines){bsw.println(new ByteBuilder(new String(nl)));}   //embed net verbatim
		final float[] sc=new float[nf]; final byte[] lb=new byte[nf];
		for(int L=0; L<NL; L++){
			for(int i=0; i<nf; i++){sc[i]=scoreRows.get(i)[L]; lb[i]=labRows.get(i)[L];}
			final float[][] kn=pav(sc, lb, knots);
			bsw.println(new ByteBuilder("#label "+labels[L]));
			final ByteBuilder h=new ByteBuilder("#lut "); h.append(kn[0].length); bsw.println(h);
			for(int j=0; j<kn[0].length; j++){
				final ByteBuilder b=new ByteBuilder(); b.append(kn[0][j], 6).tab().append(kn[1][j], 6); bsw.println(b);
			}
		}
		bsw.poisonAndWait();
		System.err.println("wrote multi-output bundle: "+outPath);
	}

	/** PAV isotonic (score->P, monotone non-decreasing), then thin to <=knots knots for a compact LUT. */
	static float[][] pav(float[] score, byte[] label, int knots){
		final int n=score.length;
		final Integer[] ord=new Integer[n];
		for(int i=0; i<n; i++){ord[i]=i;}
		java.util.Arrays.sort(ord, (a,b)->Float.compare(score[a], score[b]));
		final double[] val=new double[n], w=new double[n], ss=new double[n]; int nb=0;
		for(int k=0; k<n; k++){
			final int i=ord[k];
			val[nb]=label[i]; w[nb]=1; ss[nb]=score[i]; nb++;
			while(nb>=2 && val[nb-1]<val[nb-2]-1e-12){
				final double nw=w[nb-1]+w[nb-2];
				val[nb-2]=(val[nb-1]*w[nb-1]+val[nb-2]*w[nb-2])/nw;
				ss[nb-2]=ss[nb-1]+ss[nb-2]; w[nb-2]=nw; nb--;
			}
		}
		final float[] fx=new float[nb], fy=new float[nb];
		for(int b=0; b<nb; b++){fx[b]=(float)(ss[b]/w[b]); fy[b]=(float)val[b];}
		if(nb<=knots){return new float[][]{fx, fy};}
		//thin: keep evenly-spaced knots (monotone preserved)
		final float[] kx=new float[knots], ky=new float[knots];
		for(int j=0; j<knots; j++){
			final int idx=(int)((long)j*(nb-1)/(knots-1));
			kx[j]=fx[idx]; ky[j]=fy[idx];
		}
		return new float[][]{kx, ky};
	}
}
