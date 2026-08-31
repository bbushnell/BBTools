package prok;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import consensus.BaseGraph;
import fileIO.FileFormat;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/**
 * Training-vector generator for the ncRNA boundary-precision NN: reads a flanked
 * corpus (cutgff flank=20 output), picks the best consensus model per record via
 * the family's own 7-mer index, and emits labeled feature vectors for true and
 * shifted boundary candidates. Reuses TrnaBoundaryFeatures' generic features
 * (ANI, enrichment profile, tip fuzziness, length ratio, contig GC) -- omits
 * stemFeature (acceptor-stem palindrome is tRNA-specific).
 *
 * <p>10 features: ANI(1) + enrichmentProfile(3) + isStop(1) + tipFuzziness(3) +
 * lengthRatio(1) + contigGC(1). Same #dims TSV format as TrnaBoundaryVectorGen.
 * @author Noire
 */
public class NcrnaBoundaryVectorGen {

	public static void main(String[] args) throws IOException {
		String fastaPath=null, outPath=null, tableStartPath=null, tableStopPath=null, familyName=null;
		float meanLen=0;
		for(String a : args){
			String[] kv=a.split("=", 2);
			if(kv.length<2){continue;}
			if(kv[0].equalsIgnoreCase("fasta")){fastaPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("out")){outPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestart")){tableStartPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestop")){tableStopPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("family")){familyName=kv[1];}
			else if(kv[0].equalsIgnoreCase("meanlen")){meanLen=Float.parseFloat(kv[1]);}
		}
		if(fastaPath==null || outPath==null || tableStartPath==null || tableStopPath==null || familyName==null){
			System.err.println("Usage: family=<rnasep|srp_small|srp_large|tmrna> fasta=<flanked.fa> out=<vectors.tsv>");
			System.err.println("  tablestart=<start_table.tsv> tablestop=<stop_table.tsv> [meanlen=380]");
			System.exit(1);
		}

		familyName=CallGenes.parseNcrnaFamily(familyName);
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		if(familyName.equals("tmrna")){CallGenes.TMRNA_ENABLED=true;}
		CallGenes.loadNcrnaResources();
		NcrnaFamily fam=null;
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){
			if(f.name.equals(familyName)){fam=f; break;}
		}
		if(fam==null || fam.library==null){
			System.err.println("FATAL: family '"+familyName+"' not found or has no library.");
			System.exit(2);
		}
		if(meanLen<=0){meanLen=medianLength(fam.library);}

		final byte[][] library=fam.library;
		final BaseGraph[] models=fam.models;
		final TrnaKmerIndex index=new TrnaKmerIndex(library, 7, fam.adaptive,
			fam.adaptFloor, fam.adaptTopFrac, fam.adaptQFrac, fam.fixedMinHits);

		final TrnaNinemerTableBuilder.LoadedTable lt1=TrnaNinemerTableBuilder.loadTable(tableStartPath);
		final TrnaNinemerTableBuilder.LoadedTable lt2=TrnaNinemerTableBuilder.loadTable(tableStopPath);
		assert(lt1.type==TrnaBoundaryFeatures.BoundaryType.START);
		assert(lt2.type==TrnaBoundaryFeatures.BoundaryType.STOP);
		final TrnaBoundaryFeatures.NinemerTable startTable=lt1.table, stopTable=lt2.table;
		final int startInside=lt1.insideCount, startOutside=lt1.outsideCount;
		final int stopInside=lt2.insideCount, stopOutside=lt2.outsideCount;
		System.err.println("Family: "+familyName+" ("+library.length+" models, meanLen="+meanLen
			+", indexTopN="+fam.indexTopN+")");
		System.err.println("Tables: start inside="+startInside+" outside="+startOutside
			+", stop inside="+stopInside+" outside="+stopOutside);

		Shared.TRIM_READ_DESCRIPTION=false;
		Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final ArrayList<Read> reads=ReadInputStream.toReads(fastaPath, FileFormat.FA, -1);
		assert(!reads.isEmpty()) : "Empty flanked fasta: "+fastaPath;
		System.err.println("Loaded "+reads.size()+" records.");

		long noFlank=0, malformed=0, noModel=0, noGC=0;
		int written=0;
		final float meanLenFinal=meanLen;
		try(PrintStream out=new PrintStream(new FileOutputStream(outPath))){
			out.println("#dims\t10\t1");
			for(Read r : reads){
				final int lf=TrnaNinemerTableBuilder.parseFlankValue(r.id, "lflank=");
				final int rf=TrnaNinemerTableBuilder.parseFlankValue(r.id, "rflank=");
				if(lf<0 || rf<0){noFlank++; continue;}
				final byte[] bases=r.bases;
				final int trueStart=lf;
				final int trueStop=bases.length-rf-1;
				if(trueStart<0 || trueStop>=bases.length || trueStop<=trueStart){malformed++; continue;}
				final float contigGC=parseHeaderFloat(r.id, "contig_gc=");
				if(Float.isNaN(contigGC)){noGC++; continue;}

				final byte[] geneSeq=java.util.Arrays.copyOfRange(bases, trueStart, trueStop+1);
				//Training must use the family's inference shortlist width. A hardcoded smaller
				//value creates train/serve skew by assigning some loci a different consensus model.
				final int[] shortlist=index.shortlist(geneSeq, fam.indexTopN);
				if(shortlist.length==0){noModel++; continue;}
				int bestModel=-1; float bestId=-1;
				for(int m : shortlist){
					float id=TrnaBoundaryFeatures.aniFeature(geneSeq, library[m]);
					if(id>bestId){bestId=id; bestModel=m;}
				}
				if(bestModel<0){noModel++; continue;}
				final BaseGraph model=(models!=null && bestModel<models.length ? models[bestModel] : null);

				written+=emitVectors(out, bases, trueStart, trueStop, true, library[bestModel],
					startTable, startInside, startOutside, model, contigGC, meanLenFinal, fam.boundaryStartOffsets);
				written+=emitVectors(out, bases, trueStart, trueStop, false, library[bestModel],
					stopTable, stopInside, stopOutside, model, contigGC, meanLenFinal, fam.boundaryStopOffsets);
			}
		}
		System.err.println("Wrote "+written+" vectors. Skipped: "+noFlank+" no-flank, "
			+malformed+" malformed, "+noGC+" no-gc, "+noModel+" no-model.");
	}

	private static int emitVectors(PrintStream out, byte[] window, int trueStart, int trueStop,
			boolean varyStart, byte[] modelConsensus, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, BaseGraph model, float contigGC, float meanLen, int[] offsets){
		final TrnaBoundaryFeatures.BoundaryType type=(varyStart
			? TrnaBoundaryFeatures.BoundaryType.START : TrnaBoundaryFeatures.BoundaryType.STOP);
		final int trueBoundaryPos=(varyStart ? trueStart : trueStop);
		int n=0;
		for(int offset : offsets){
			int s=(varyStart ? trueStart+offset : trueStart);
			int e=(varyStart ? trueStop : trueStop+offset);
			if(s<0 || e>=window.length || e-s<15){continue;}
			int label=(offset==0 ? 1 : 0);
			int boundaryPos=(varyStart ? s : e);
			final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
			float ani=TrnaBoundaryFeatures.aniFeature(candSeq, modelConsensus);
			float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, trueBoundaryPos,
				type, insideCount, outsideCount, table);
			float isStop=(varyStart ? 0f : 1f);
			float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, varyStart);
			float lengthRatio=(e-s+1)/meanLen;
			out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",
				ani, prof[0], prof[1], prof[2], isStop, fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC, label);
			n++;
		}
		return n;
	}

	private static float medianLength(byte[][] library){
		int[] lens=new int[library.length];
		for(int i=0; i<library.length; i++){lens[i]=library[i].length;}
		java.util.Arrays.sort(lens);
		return lens[lens.length/2];
	}

	private static float parseHeaderFloat(String header, String key){
		final int idx=header.indexOf(key);
		if(idx<0){return Float.NaN;}
		final int start=idx+key.length();
		int end=start;
		while(end<header.length() && (Character.isDigit(header.charAt(end)) || header.charAt(end)=='.')){end++;}
		if(end==start){return Float.NaN;}
		try{return Float.parseFloat(header.substring(start, end));}
		catch(NumberFormatException e){return Float.NaN;}
	}
}
