package var2;

import java.io.PrintStream;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Copies neural-network training labels from a labeled "master" VCF (the output of
 * {@link LabelVCF}, computed on the deepest-coverage data) onto a different VCF of the SAME
 * variants called at a DIFFERENT coverage.  This is the mechanism behind depth-independent
 * labeling: the label is computed once on the deepest data and transferred to every lower-
 * coverage call set, so the same variant carries the same target at 10x/20x/40x/80x while its
 * feature vector differs by depth.
 *
 * Both the master and the input VCF must be CallVariants output on the SAME reference, so a
 * variant's identity key CHROM:POS:REF:ALT is identical between them (only the support metrics
 * differ with coverage).  No normalization or scaffold map is required.
 *
 * For each input variant: if its key is present in the master, the master's ";NNL=...;NNC=..."
 * suffix is spliced onto the input line's INFO column.  Input variants absent from the master
 * (low-coverage-only artifacts that never reached count>=2 at full depth) are emitted UNLABELED;
 * the downstream vector generator drops any line lacking NNL.
 *
 * Usage: java var2.TransferLabels master=labeled_80x.vcf.gz in=allvariants_40x.vcf.gz out=labeled_40x.vcf.gz
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class TransferLabels {

	public static void main(String[] args){
		String masterFile=null, inFile=null, outFile=null;

		for(String arg : args){
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("master") || a.equals("labels")){masterFile=b;}
			else if(a.equals("in")){inFile=b;}
			else if(a.equals("out")){outFile=b;}
		}

		assert(masterFile!=null) : "Missing master= (labeled master VCF from LabelVCF)";
		assert(inFile!=null) : "Missing in= (a CallVariants VCF to receive labels)";
		assert(outFile!=null) : "Missing out= (labeled output VCF)";

		PrintStream out=System.err;

		// 1) Load the master labels: key CHROM:POS:REF:ALT -> the ";NNL=...;NNC=..." INFO suffix.
		out.println("Loading master labels: "+masterFile);
		HashMap<String, String> labels=new HashMap<String, String>();
		ByteFile mbf=ByteFile.makeByteFile(FileFormat.testInput(masterFile, FileFormat.TXT, null, true, false));
		long masterRows=0, masterNoLabel=0;
		for(byte[] line=mbf.nextLine(); line!=null; line=mbf.nextLine()){
			if(line.length>0 && line[0]=='#'){continue;}
			String[] f=new String(line).split("\t");
			if(f.length<8){continue;}
			int idx=f[7].indexOf(";NNL=");
			if(idx<0){masterNoLabel++; continue;}
			String key=f[0]+":"+f[1]+":"+f[3]+":"+f[4];
			labels.put(key, f[7].substring(idx)); // ";NNL=x;NNC=y" (to end of INFO)
			masterRows++;
		}
		mbf.close();
		out.println("Master labels loaded: "+masterRows+" (master rows without NNL: "+masterNoLabel+")");

		// 2) Stream the input VCF, splice matching labels into INFO, emit.
		out.println("Transferring onto: "+inFile);
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(inFile, FileFormat.TXT, null, true, false));
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		long total=0, transferred=0, unmatched=0;
		boolean headerInjected=false;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]=='#'){
				if(!headerInjected && Tools.startsWith(line, "#CHROM")){
					bb.clear();
					bb.append("##INFO=<ID=NNL,Number=1,Type=Float,Description=\"UMP45 NN training target label (transferred)\">").nl();
					bb.append("##INFO=<ID=NNC,Number=1,Type=Integer,Description=\"UMP45 label category: 1=GIAB-HC,2=GIAB-LC,3=nonGIAB-HCregion,4=nonGIAB-LCregion\">").nl();
					bsw.print(bb);
					headerInjected=true;
				}
				bb.clear(); bb.append(line).nl(); bsw.print(bb);
				continue;
			}

			String[] f=new String(line).split("\t");
			if(f.length<8){bb.clear(); bb.append(line).nl(); bsw.print(bb); continue;}
			total++;
			String key=f[0]+":"+f[1]+":"+f[3]+":"+f[4];
			String suffix=labels.get(key);

			bb.clear();
			for(int i=0; i<f.length; i++){
				if(i>0){bb.tab();}
				bb.append(f[i]);
				if(i==7 && suffix!=null){bb.append(suffix);}
			}
			bb.nl();
			bsw.print(bb);

			if(suffix!=null){transferred++;}else{unmatched++;}
		}
		bf.close();
		bsw.poisonAndWait();

		out.println("Input variants:        "+total);
		out.println("  labels transferred:  "+transferred);
		out.println("  unmatched (no label):"+unmatched);
		out.println("Output: "+outFile);
	}
}
