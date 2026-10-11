package var2;

import java.io.PrintStream;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Generates neural-network training vectors from a label-annotated VCF. The
 * label usually comes from {@link LabelVCF} and is stored in the NNL INFO key.
 * This is the depth-independent half of the pipeline: the NN target rides with
 * the variant, while the feature vector is recomputed from this VCF's depth.
 * The same variant therefore gets the same target but a different vector at
 * 10x, 40x and 160x.
 *
 * The VCF is streamed one line at a time. Each variant is parsed with
 * {@link VcfToVar#fromVCF}, its label is read from the semicolon-delimited NNL
 * key, and the {@link VectorUMP45} feature vector is emitted followed by the
 * target. Variants whose label falls in the deadzone (0.5 +/- deadzone, default
 * the 0.4-0.6 band) are excluded so ambiguous examples never reach the net.
 * Variant rows are streamed rather than retained, while scaffold metadata and
 * optional reference bases remain resident. The vector math is reused verbatim
 * from VectorUMP45; only the label source differs from VcfToTrainingVectors
 * (read, not recomputed).
 *
 * Dataset-level statistics (pairing rate, quality/mapq/readlen averages,
 * ploidy) are read from exact ## header keys matching VcfLoader output. Without
 * ref=, the VCF ##contig header supplies only names and lengths; sequence-based
 * HMP and N-run CED features then fall back to no-bases values.
 *
 * Usage: java var2.VectorsFromLabeledVCF in=labeled.vcf.gz out=vectors.tsv [ref=ref.fa.gz]
 *        [deadzone=0.1] [ploidy=N] [platform=illumina] [includescore=f]
 *
 * @author UMP45
 * @contributor Sucrose
 * @date June 7, 2026
 */
public class VectorsFromLabeledVCF {

	public static void main(String[] args){
		String inFile=null, refFile=null, outFile=null;
		double deadzone=0.1;        // half-width of the excluded band around 0.5 (0.1 => skip 0.4-0.6)
		int ploidyOverride=-1;      // if >0, overrides the ##ploidy from the header
		HashSet<String> allowedChroms=null; // if set, only these CHROM values are vectorized

		for(String arg : args){
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in")){inFile=b;}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("deadzone")){deadzone=Double.parseDouble(b);}
			else if(a.equals("ploidy")){ploidyOverride=Integer.parseInt(b);}
			else if(a.equals("includescore")){VectorUMP45.includeScore=Parse.parseBoolean(b);}
			else if(a.equals("scaffolds") || a.equals("chroms")){
				allowedChroms=new HashSet<String>();
				for(String c : b.split(",")){allowedChroms.add(c.trim());}
			}
			else if(a.equals("platform")){
				if(b.equalsIgnoreCase("illumina")){VectorUMP45.platform=VectorUMP45.PLATFORM_ILLUMINA;}
				else if(b.equalsIgnoreCase("pacbio")){VectorUMP45.platform=VectorUMP45.PLATFORM_PACBIO;}
				else if(b.equalsIgnoreCase("nanopore") || b.equalsIgnoreCase("ont")){VectorUMP45.platform=VectorUMP45.PLATFORM_NANOPORE;}
				else if(b.equalsIgnoreCase("roche") || b.equalsIgnoreCase("sbx")){VectorUMP45.platform=VectorUMP45.PLATFORM_ROCHE;}
				else{VectorUMP45.platform=Integer.parseInt(b);}
			}
		}

		assert(inFile!=null) : "Missing in= (labeled VCF)";
		assert(outFile!=null) : "Missing out= (vector TSV)";
		assert(deadzone>=0 && deadzone<0.5) : "deadzone must be in [0,0.5): "+deadzone;

		PrintStream out=System.err;

		// Scaffold map for VcfToVar.fromVCF / makeVector.
		//TODO: Possible bug [var2/VectorsFromLabeledVCF#002] (MEDIUM, feature drift) - without ref=,
		//VectorUMP45 computes HMP as 0 and CED from scaffold length only; VcfToVar.fromVCF does not
		//parse the HMP/CED INFO values already present in the VCF. Either require ref= or cache the
		//VCF-provided HMP/CED fields for post-hoc vector generation.
		ScafMap scafMap;
		if(refFile!=null){
			out.println("Loading reference: "+refFile);
			scafMap=ScafMap.loadReference(refFile, true);
		}else{
			out.println("Loading scaffold map from VCF header: "+inFile);
			scafMap=ScafMap.loadVcfHeader(inFile);
		}

		// Dataset-level stats (overwritten from ## header lines as we stream).
		double pairingRate=0, totalQualityAvg=0, totalMapqAvg=0, readLengthAvg=0;
		int ploidy=2;

		final double dzLo=0.5-deadzone, dzHi=0.5+deadzone;

		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(inFile, FileFormat.TXT, null, true, false));
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		bb.append("#dims\t").append(VectorUMP45.DIMS).append('\t').append(1).nl();
		bsw.print(bb);

		long written=0, excluded=0, positive=0, negative=0, nolabel=0, filtered=0;

		Throwable pending=null;
		boolean errorState=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0 && line[0]=='#'){
					// Capture dataset stats from ## header lines (key=value).
					String h=new String(line);
					int eq=h.indexOf('=');
					if(eq>0){
						String k=h.substring(0, eq), val=h.substring(eq+1).trim();
						if(k.equalsIgnoreCase("##properPairRate")){pairingRate=Double.parseDouble(val);}
						else if(k.equalsIgnoreCase("##totalQualityAvg")){totalQualityAvg=Double.parseDouble(val);}
						else if(k.equalsIgnoreCase("##mapqAvg")){totalMapqAvg=Double.parseDouble(val);}
						else if(k.equalsIgnoreCase("##readLengthAvg")){readLengthAvg=Double.parseDouble(val);}
						else if(k.equalsIgnoreCase("##ploidy")){ploidy=Integer.parseInt(val);}
					}
					continue;
				}

				String s=new String(line);
				if(allowedChroms!=null){
					int t=s.indexOf('\t');
					String chrom=(t>0 ? s.substring(0, t) : s);
					if(!allowedChroms.contains(chrom)){filtered++; continue;}
				}
				int infoStart=-1, infoEnd=s.length(), tabs=0;
				for(int i=0; i<s.length() && tabs<8; i++){
					if(s.charAt(i)=='\t'){
						tabs++;
						if(tabs==7){infoStart=i+1;}
						else if(tabs==8){infoEnd=i;}
					}
				}
				if(infoStart<0){nolabel++; continue;}
				int idx=s.startsWith("NNL=", infoStart) ? infoStart : s.indexOf(";NNL=", infoStart);
				if(idx<0){nolabel++; continue;}
				if(idx>=infoEnd){nolabel++; continue;}
				if(s.charAt(idx)==';'){idx++;}
				int lstart=idx+4;
				int lend=s.indexOf(';', lstart);
				if(lend<0 || lend>infoEnd){lend=infoEnd;}
				double label=Double.parseDouble(s.substring(lstart, lend));

				// Deadzone exclusion (ambiguous labels never train).
				if(label>dzLo && label<dzHi){excluded++; continue;}

				// Parse the variant and compute its feature vector.
				Var v;
				try{
					v=VcfToVar.fromVCF(line, scafMap, true, true);
				}catch(Exception e){
					nolabel++; continue;
				}
				if(v==null){nolabel++; continue;}
				v.calcCoverage(scafMap);

				final int useploidy=(ploidyOverride>0 ? ploidyOverride : ploidy);
				float[] vec=VectorUMP45.makeVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, useploidy, scafMap);

				bb.clear();
				for(int i=0; i<vec.length; i++){
					if(i>0){bb.tab();}
					bb.append(vec[i], 6);
				}
				bb.tab().append(label, 6).nl();
				bsw.print(bb);

				written++;
				if(label>=0.5){positive++;}else{negative++;}
			}
		}catch(Throwable t){
			pending=t;
		}
		// Close the reader and writer on normal and exceptional exits.
		try{
			errorState=bf.close();
		}catch(Throwable t){
			pending=suppressOrReturn(pending, t);
		}
		try{
			errorState|=bsw.poisonAndWait();
		}catch(Throwable t){
			pending=suppressOrReturn(pending, t);
		}
		if(errorState){
			pending=suppressOrReturn(pending, new RuntimeException(VectorsFromLabeledVCF.class.getName()+
					" terminated in an error state; the output may be corrupt."));
		}
		if(pending!=null){rethrow(pending);}

		out.println("Vectors written: "+written+" (positive: "+positive+", negative: "+negative+")");
		out.println("Excluded (deadzone "+dzLo+"-"+dzHi+"): "+excluded);
		out.println("Filtered (scaffold not in include list): "+filtered);
		out.println("Skipped (no label / parse fail): "+nolabel);
		out.println("Output: "+outFile);
	}

	/** Returns pending with the new cleanup failure suppressed, preserving the primary cause. */
	private static Throwable suppressOrReturn(final Throwable pending, final Throwable cleanup){
		if(pending==null){return cleanup;}
		pending.addSuppressed(cleanup);
		return pending;
	}

	/** Rethrows without wrapping unchecked exceptions or errors. */
	private static void rethrow(final Throwable t){
		if(t instanceof RuntimeException){throw (RuntimeException)t;}
		if(t instanceof Error){throw (Error)t;}
		throw new RuntimeException(t);
	}
}
