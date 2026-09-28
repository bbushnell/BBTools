package prok;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import consensus.BaseGraph;
import fileIO.FileFormat;
import idaligner.QuantumAligner;
import map.LongHashSet;
import parse.Parse;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/**
 * Training-vector generator for the ncRNA boundary-precision NN: reads a flanked
 * corpus (cutgff flank=20 output), picks the best consensus model per record via
 * the family's configured mapping-k index, and emits labeled feature vectors for true and
 * shifted boundary candidates. Reuses TrnaBoundaryFeatures' generic features
 * (ANI, enrichment profile, tip fuzziness, length ratio, contig GC) -- omits
 * stemFeature (acceptor-stem palindrome is tRNA-specific).
 *
 * <p>V3's 10 features are enrichmentProfile(3), isStop(1), lengthRatio(1),
 * contigGC(1), accepted-locus Quantum identity(1), and three explicit zero spare
 * slots. V1/V2 retain their frozen layouts for reproducibility.
 * @author Noire
 */
public class NcrnaBoundaryVectorGen {

	public static void main(String[] args) throws IOException {
		String fastaPath=null, outPath=null, tableStartPath=null, tableStopPath=null, familyTableStartPath=null, familyTableStopPath=null, familyName=null;
		String statsPath=null, boundaryFuzzPath=null, fixedModelName=null, trainerOutPath=null;
		int boundaryFeatureVersion=NcrnaBoundaryScorer.FEATURES_V1;
		String r58Kmers=null, r58Consensus=null, r58Models=null;
		String lsuKmers=null, lsuConsensus=null, lsuModels=null;
		String externalName=null, externalConsensus=null, externalModels=null, externalSeeds=null;
		int externalK=-1, externalSeedMinHits=-1, externalMinLen=-1, externalWindowPad=-1;
		int externalIndexK=-1, externalFixedMinHits=-1, externalTopN=-1;
		float externalIdPass=Float.NaN, externalIdBorderline=Float.NaN;
		Boolean externalAdaptive=null, externalRankedFallback=null, externalStrictIndexCutoff=null;
		int[] externalStartOffsets=null, externalStopOffsets=null;
		boolean stopOnly=false;
		boolean leaveOneOut=true;
		boolean countOnly=false;
		boolean modelColumn=false;
		float meanLen=0;
		for(String a : args){
			String[] kv=a.split("=", 2);
			if(kv.length<2 || kv[0].isEmpty() || kv[1].isEmpty()){
				throw new IllegalArgumentException("Expected nonempty flag=value for vector generation: "+a);
			}
			if(kv[0].equalsIgnoreCase("fasta")){fastaPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("out")){outPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestart")){tableStartPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestop")){tableStopPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("familytablestart")){familyTableStartPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("familytablestop")){familyTableStopPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("statsout")){statsPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("family")){familyName=kv[1];}
			else if(kv[0].equalsIgnoreCase("boundaryfeatures")){boundaryFeatureVersion=NcrnaBoundaryScorer.parseFeatureVersion(kv[1]);}
			else if(kv[0].equalsIgnoreCase("boundaryfuzz") || kv[0].equalsIgnoreCase("boundaryfuzzconstants")){boundaryFuzzPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("fixedmodel")){fixedModelName=kv[1];}
			else if(kv[0].equalsIgnoreCase("trainerout")){trainerOutPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("modelcolumn")){modelColumn=Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("meanlen")){meanLen=CallGenes.parseFiniteFloat("meanlen",kv[1]);}
			else if(kv[0].equalsIgnoreCase("r58kmers")){r58Kmers=kv[1];}
			else if(kv[0].equalsIgnoreCase("r58consensus")){r58Consensus=kv[1];}
			else if(kv[0].equalsIgnoreCase("r58models")){r58Models=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsukmers")){lsuKmers=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsuconsensus")){lsuConsensus=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsumodels")){lsuModels=kv[1];}
			else if(kv[0].equalsIgnoreCase("externalname")){externalName=kv[1];}
			else if(kv[0].equalsIgnoreCase("externalconsensus")){externalConsensus=kv[1];}
			else if(kv[0].equalsIgnoreCase("externalmodels")){externalModels=kv[1];}
			else if(kv[0].equalsIgnoreCase("externalseeds")){externalSeeds=kv[1];}
			else if(kv[0].equalsIgnoreCase("externalk")){externalK=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalseedminhits")){externalSeedMinHits=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalminlen")){externalMinLen=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalwindowpad")){externalWindowPad=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalindexk")){externalIndexK=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalfixedminhits")){externalFixedMinHits=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externaltopn")){externalTopN=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalidpass")){externalIdPass=Float.parseFloat(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalidborderline")){externalIdBorderline=Float.parseFloat(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externaladaptive")){externalAdaptive=Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalrankedfallback")){externalRankedFallback=Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalstrictindexcutoff")){externalStrictIndexCutoff=Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("externalstartoffsets")){externalStartOffsets=parseOffsets(kv[1], "start");}
			else if(kv[0].equalsIgnoreCase("externalstopoffsets")){externalStopOffsets=parseOffsets(kv[1], "stop");}
			else if(kv[0].equalsIgnoreCase("stoponly") || kv[0].equalsIgnoreCase("3primeonly")){
				stopOnly=Parse.parseBoolean(kv[1]);
			}else if(kv[0].equalsIgnoreCase("leaveoneout") || kv[0].equalsIgnoreCase("loo")){
				leaveOneOut=Parse.parseBoolean(kv[1]);
			}else if(kv[0].equalsIgnoreCase("countonly")){
				countOnly=Parse.parseBoolean(kv[1]);
			}else{throw new IllegalArgumentException("Unknown vector option: "+kv[0]);}
		}
		if(fastaPath==null || outPath==null || tableStartPath==null || tableStopPath==null || familyName==null){
			System.err.println("Usage: family=<rnasep|srp_small|srp_large|tmrna> fasta=<flanked.fa> out=<vectors.tsv>");
			System.err.println("  tablestart=<start_table.tsv> tablestop=<stop_table.tsv> [meanlen=380] [stoponly=f]");
			System.exit(1);
		}
		if(countOnly && statsPath==null){
			throw new IllegalArgumentException("countonly=t requires statsout= so the source-derived accounting is durable");
		}
		if(modelColumn && (fixedModelName==null || trainerOutPath==null)){
			throw new IllegalArgumentException("modelcolumn=t requires fixedmodel= and trainerout=");
		}
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V3 && (familyTableStartPath==null || familyTableStopPath==null)){
			throw new IllegalArgumentException("boundaryfeatures=v3 requires familytablestart= and familytablestop= for the 50/50 family/consensus profile blend");
		}

		final boolean anyExternal=(externalName!=null || externalConsensus!=null || externalModels!=null || externalSeeds!=null ||
			externalK>=0 || externalSeedMinHits>=0 || externalMinLen>=0 || externalWindowPad>=0 ||
			externalIndexK>=0 || externalFixedMinHits>=0 || externalTopN>=0 ||
			!Float.isNaN(externalIdPass) || !Float.isNaN(externalIdBorderline) || externalAdaptive!=null ||
			externalRankedFallback!=null || externalStrictIndexCutoff!=null ||
			externalStartOffsets!=null || externalStopOffsets!=null);
		final boolean externalMode=familyName.equalsIgnoreCase("external");
		if(anyExternal && !externalMode){
			throw new IllegalArgumentException("Explicit external-family parameters cannot be combined with named built-in family='"
				+familyName+"'; use family=external so built-in resource loading cannot be silently overridden.");
		}
		NcrnaFamily fam=null;
		if(externalMode){
			fam=loadExternalFamily(externalName, externalConsensus, externalModels, externalSeeds, externalK, externalSeedMinHits,
				externalMinLen, externalWindowPad, externalIndexK, externalFixedMinHits, externalTopN,
				externalIdPass, externalIdBorderline, externalAdaptive, externalRankedFallback,
				externalStrictIndexCutoff, externalStartOffsets, externalStopOffsets);
			familyName=fam.name;
		}else{
			familyName=CallGenes.parseNcrnaFamily(familyName);
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			if(familyName.equals("tmrna")){CallGenes.TMRNA_ENABLED=true;}
			if(familyName.equals("r58") || familyName.equals("lsu")){
				CallGenes.setR58LsuEnabled(true);
				CallGenes.R58_KMERS_OVERRIDE=r58Kmers;
				CallGenes.R58_CONSENSUS_OVERRIDE=r58Consensus;
				CallGenes.R58_MODELS_OVERRIDE=r58Models;
				CallGenes.LSU_KMERS_OVERRIDE=lsuKmers;
				CallGenes.LSU_CONSENSUS_OVERRIDE=lsuConsensus;
				CallGenes.LSU_MODELS_OVERRIDE=lsuModels;
			}
			CallGenes.loadNcrnaResources();
			for(NcrnaFamily f : GeneCaller.ncrnaFamilies){
				if(f.name.equals(familyName)){fam=f; break;}
			}
		}
		if(fam==null || fam.library==null){
			System.err.println("FATAL: family '"+familyName+"' not found or has no library.");
			System.exit(2);
		}
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2 && boundaryFuzzPath==null){
			throw new IllegalArgumentException("boundaryfeatures=v2 requires boundaryfuzz=<model/end constants TSV>");
		}
		if(meanLen<=0){meanLen=medianLength(fam.library);}

		final byte[][] library=fam.library;
		final BaseGraph[] models=fam.models;
		int fixedModel=-1;
		if(fixedModelName!=null){for(int i=0;i<fam.modelNames.length;i++){if(firstToken(fam.modelNames[i]).equals(fixedModelName)){fixedModel=i;break;}}if(fixedModel<0){throw new IllegalArgumentException("fixedmodel is absent from consensus library: "+fixedModelName);}}
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2 && (models==null || models.length!=library.length)){
			throw new IllegalArgumentException("boundaryfeatures=v2 requires one HBM model per consensus");
		}
		final NcrnaBoundaryFuzzConstants boundaryFuzz=(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2
			? NcrnaBoundaryFuzzConstants.load(boundaryFuzzPath,fam.modelNames) : null);
		//Ranking uses the same family setting as inference; endpoint-table K is separate.
		final int effectiveIndexK=fam.indexK;
		final TrnaKmerIndex index=new TrnaKmerIndex(library, effectiveIndexK, fam.adaptive,
			fam.adaptFloor, fam.adaptTopFrac, fam.adaptQFrac, fam.fixedMinHits);

		final TrnaNinemerTableBuilder.LoadedTable lt1=TrnaNinemerTableBuilder.loadTable(tableStartPath);
		final TrnaNinemerTableBuilder.LoadedTable lt2=TrnaNinemerTableBuilder.loadTable(tableStopPath);
		assert(lt1.type==TrnaBoundaryFeatures.BoundaryType.START);
		assert(lt2.type==TrnaBoundaryFeatures.BoundaryType.STOP);
		final TrnaBoundaryFeatures.NinemerTable startTable=lt1.table, stopTable=lt2.table;
		final int startInside=lt1.insideCount, startOutside=lt1.outsideCount;
		final int stopInside=lt2.insideCount, stopOutside=lt2.outsideCount;
		final TrnaBoundaryFeatures.NinemerTable familyStartTable,familyStopTable;
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V3){
			final TrnaNinemerTableBuilder.LoadedTable familyStart=TrnaNinemerTableBuilder.loadTable(familyTableStartPath),familyStop=TrnaNinemerTableBuilder.loadTable(familyTableStopPath);
			validateFamilyConsensusTable(lt1,familyStart,TrnaBoundaryFeatures.BoundaryType.START,"start");
			validateFamilyConsensusTable(lt2,familyStop,TrnaBoundaryFeatures.BoundaryType.STOP,"stop");
			familyStartTable=familyStart.table;familyStopTable=familyStop.table;
		}else{familyStartTable=null;familyStopTable=null;}
		System.err.println("Family: "+familyName+" ("+library.length+" models, meanLen="+meanLen
			+", indexK="+effectiveIndexK+", indexTopN="+fam.indexTopN+", boundaryFeatures=v"+boundaryFeatureVersion+")");
		System.err.println("Tables: start inside="+startInside+" outside="+startOutside
			+", stop inside="+stopInside+" outside="+stopOutside+", leaveOneOut="+leaveOneOut);

		Shared.TRIM_READ_DESCRIPTION=false;
		Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final ArrayList<Read> reads=ReadInputStream.toReads(fastaPath, FileFormat.FA, -1);
		assert(!reads.isEmpty()) : "Empty flanked fasta: "+fastaPath;
		System.err.println("Loaded "+reads.size()+" records.");

		long noFlank=0, malformed=0, noModel=0, noGC=0;
		int written=0;
		final VectorStats vectorStats=new VectorStats();
		vectorStats.records=reads.size();
		final float meanLenFinal=meanLen;
		final boolean quantumFeatures=(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2
			|| boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V3);
		final int[] quantumPos=(quantumFeatures ? new int[4] : null);
		final PrintStream out;
		final PrintStream trainerOut;
		try{out=new PrintStream(new FileOutputStream(outPath));trainerOut=(trainerOutPath==null ? null : new PrintStream(new FileOutputStream(trainerOutPath)));}catch(IOException e){throw new RuntimeException("Could not open vector output",e);}
		try{
			if(!countOnly){out.println(modelColumn ? "#dims\t10\t1\tmodel" : "#dims\t10\t1");if(trainerOut!=null){trainerOut.println("#dims\t10\t1");}}
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
				final int[] shortlist=(fixedModel>=0 ? new int[]{fixedModel} : index.shortlist(geneSeq, fam.indexTopN));
				if(shortlist.length==0){noModel++; continue;}
				if(fixedModel>=0){final String tagged=NcrnaConsensusPartitioner.recordedModel(r.id);if(tagged==null||!tagged.equals(fixedModelName)){throw new IllegalArgumentException("fixedmodel/header assignment mismatch: fixed="+fixedModelName+" header="+r.id);}}
				int bestModel=-1; float bestId=-1;
				for(int m : shortlist){
					final float id=(quantumFeatures
						? QuantumAligner.alignStatic(library[m],bases,quantumPos)
						: TrnaBoundaryFeatures.aniFeature(geneSeq, library[m]));
					if(id>bestId){bestId=id; bestModel=m;}
					if(quantumFeatures && externalMode
							&& Boolean.TRUE.equals(externalRankedFallback) && id>=fam.idPass){break;}
				}
				if(bestModel<0){noModel++; continue;}
				vectorStats.processedLoci++;
				final BaseGraph model=(models!=null && bestModel<models.length ? models[bestModel] : null);

				if(!stopOnly){
					written+=emitVectors(out, bases, trueStart, trueStop, true, library[bestModel],
						startTable, startInside, startOutside, model, contigGC, meanLenFinal,
						fam.boundaryStartOffsets, leaveOneOut, vectorStats, !countOnly,
						boundaryFeatureVersion,bestId,(boundaryFuzz==null ? null : boundaryFuzz.values(bestModel,true)),familyStartTable,trainerOut,firstToken(fam.modelNames[bestModel]),modelColumn);
				}
				written+=emitVectors(out, bases, trueStart, trueStop, false, library[bestModel],
					stopTable, stopInside, stopOutside, model, contigGC, meanLenFinal,
					fam.boundaryStopOffsets, leaveOneOut, vectorStats, !countOnly,
					boundaryFeatureVersion,bestId,(boundaryFuzz==null ? null : boundaryFuzz.values(bestModel,false)),familyStopTable,trainerOut,firstToken(fam.modelNames[bestModel]),modelColumn);
			}
			if(out.checkError()||(trainerOut!=null&&trainerOut.checkError())){throw new IOException("Vector write error");}
		}finally{out.close();if(trainerOut!=null){trainerOut.close();}}
		vectorStats.noFlank=noFlank;
		vectorStats.malformed=malformed;
		vectorStats.noGC=noGC;
		vectorStats.noModel=noModel;
		if(written!=vectorStats.emitted){
			throw new IllegalStateException("Vector return/accounting mismatch: written="+written
				+" emitted="+vectorStats.emitted);
		}
		if(statsPath!=null){
			vectorStats.validate(stopOnly, fam.boundaryStartOffsets, fam.boundaryStopOffsets);
			vectorStats.write(statsPath, stopOnly, fam.boundaryStartOffsets, fam.boundaryStopOffsets);
		}
		System.err.println("Wrote "+written+" vectors. Skipped: "+noFlank+" no-flank, "
			+malformed+" malformed, "+noGC+" no-gc, "+noModel+" no-model.");
		if(statsPath!=null){
			System.err.println("VectorAccounting: candidates="+vectorStats.candidates+" emitted="+vectorStats.emitted
				+" skipped="+vectorStats.skipped()+" zeroEmitted="+vectorStats.zeroEmitted
				+" skipStartBeforeWindow="+vectorStats.skipStartBeforeWindow
				+" skipStopAfterWindow="+vectorStats.skipStopAfterWindow
				+" skipShort="+vectorStats.skipShort+" stats="+statsPath);
		}
	}

	static void validateFamilyConsensusTable(TrnaNinemerTableBuilder.LoadedTable consensus,
			TrnaNinemerTableBuilder.LoadedTable family,
			TrnaBoundaryFeatures.BoundaryType expectedType, String end){
		if(consensus.type!=expectedType||family.type!=expectedType){
			throw new IllegalArgumentException("Family/consensus "+end+" profile table boundary type mismatch");
		}
		if(family.k!=consensus.k){
			throw new IllegalArgumentException("Family/consensus "+end+" profile table k mismatch: family="+family.k+", consensus="+consensus.k);
		}
		if(family.insideCount!=consensus.insideCount||family.outsideCount!=consensus.outsideCount){
			throw new IllegalArgumentException("Family/consensus "+end+" profile table geometry mismatch");
		}
	}

	/** Builds one explicitly supplied development family without touching CallGenes' built-in registry. */
	static NcrnaFamily loadExternalFamily(String name, String consensusPath, String modelPath, String seedPath,
			int kLong, int seedMinHits, int minLen, int windowPad, int indexK, int fixedMinHits, int topN,
			float idPass, float idBorderline, Boolean adaptive, Boolean rankedFallback,
			Boolean strictIndexCutoff, int[] startOffsets, int[] stopOffsets){
		final StringBuilder missing=new StringBuilder();
		if(name==null || name.isEmpty()){appendMissing(missing, "externalname");}
		if(consensusPath==null || !new File(consensusPath).isFile()){appendMissing(missing, "externalconsensus");}
		if(modelPath==null || !new File(modelPath).isFile()){appendMissing(missing, "externalmodels");}
		if(seedPath==null || !new File(seedPath).isFile()){appendMissing(missing, "externalseeds");}
		if(kLong<1 || kLong>31){appendMissing(missing, "externalk");}
		if(seedMinHits<1){appendMissing(missing, "externalseedminhits");}
		if(minLen<1){appendMissing(missing, "externalminlen");}
		if(windowPad<0){appendMissing(missing, "externalwindowpad");}
		if(indexK<1 || indexK>15){appendMissing(missing, "externalindexk");}
		if(fixedMinHits<0){appendMissing(missing, "externalfixedminhits");}
		if(topN<1){appendMissing(missing, "externaltopn");}
		if(!Float.isFinite(idPass) || idPass<=0f || idPass>1f){appendMissing(missing, "externalidpass");}
		if(!Float.isFinite(idBorderline) || idBorderline<=0f || idBorderline>idPass){appendMissing(missing, "externalidborderline");}
		if(adaptive==null || adaptive.booleanValue()){appendMissing(missing, "externaladaptive=f");}
		if(rankedFallback==null){appendMissing(missing, "externalrankedfallback");}
		if(strictIndexCutoff==null){appendMissing(missing, "externalstrictindexcutoff");}
		if(startOffsets==null){appendMissing(missing, "externalstartoffsets");}
		if(stopOffsets==null){appendMissing(missing, "externalstopoffsets");}
		if(missing.length()>0){throw new IllegalArgumentException("family=external requires valid explicit resources/config: "+missing);}
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensusPath);
		if(library==null || library.length<1){throw new IllegalArgumentException("External consensus library is empty: "+consensusPath);}
		final String[] modelNames=TrnaConsensusBuilder.lastLoadedNames;
		final BaseGraph[] models=TrnaConsensusBuilder.loadModels(modelPath);
		if(models==null || models.length!=library.length || modelNames==null || modelNames.length!=library.length){
			throw new IllegalArgumentException("External consensus/HBM count mismatch: consensus="+library.length
				+" models="+(models==null ? -1 : models.length)+" names="+(modelNames==null ? -1 : modelNames.length));
		}
		for(int i=0; i<library.length; i++){
			if(!firstToken(modelNames[i]).equals(firstToken(models[i].name))){
				throw new IllegalArgumentException("External consensus/HBM name mismatch at index "+i+": fasta='"
					+modelNames[i]+"' hbm='"+models[i].name+"'");
			}
		}
		final LongHashSet seeds=ProkObject.loadLongKmers(seedPath, kLong);
		if(seeds==null || seeds.size()<1){throw new IllegalArgumentException("External seed set is empty: "+seedPath);}
		System.err.println("External family: "+name+" consensus="+consensusPath+" models="+modelPath+" seeds="+seedPath
			+" k="+kLong+" seedMinHits="+seedMinHits+" indexK="+indexK+" fixedMinHits="+fixedMinHits
			+" topN="+topN+" idPass="+idPass+" idBorderline="+idBorderline+" minLen="+minLen
			+" windowPad="+windowPad+" rankedFallback="+rankedFallback+" strictIndexCutoff="+strictIndexCutoff);
		return new NcrnaFamily(name, library, models, modelNames, seeds, kLong, minLen, windowPad,
			indexK, topN, false, 0f, 0f, 0f, fixedMinHits, 0f, 1f, idPass, idBorderline,
			0f, 0.85f, null, null, null, null, -1, -1, -1, -1, 0f,
			startOffsets, stopOffsets, 0f, 0f);
	}

	private static void appendMissing(StringBuilder sb, String key){
		if(sb.length()>0){sb.append(", ");}
		sb.append(key);
	}

	private static String firstToken(String value){
		if(value==null){return "";}
		int end=0; while(end<value.length() && !Character.isWhitespace(value.charAt(end))){end++;}
		return value.substring(0, end);
	}

	static int[] parseOffsets(String value, String label){
		final String[] split=value.split(",");
		if(split.length<1){throw new IllegalArgumentException("Empty external "+label+" offsets");}
		final int[] offsets=new int[split.length];
		for(int i=0; i<split.length; i++){offsets[i]=Integer.parseInt(split[i].trim());}
		return offsets;
	}

	static int emitVectors(PrintStream out, byte[] window, int trueStart, int trueStop,
			boolean varyStart, byte[] modelConsensus, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, BaseGraph model, float contigGC, float meanLen,
			int[] offsets, boolean leaveOneOut, VectorStats stats, boolean writeVectors){
		return emitVectors(out,window,trueStart,trueStop,varyStart,modelConsensus,table,insideCount,outsideCount,
			model,contigGC,meanLen,offsets,leaveOneOut,stats,writeVectors,
			NcrnaBoundaryScorer.FEATURES_V1,Float.NaN,null);
	}

	static int emitVectors(PrintStream out, byte[] window, int trueStart, int trueStop,
			boolean varyStart, byte[] modelConsensus, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, BaseGraph model, float contigGC, float meanLen,
			int[] offsets, boolean leaveOneOut, VectorStats stats, boolean writeVectors,
			int featureVersion, float locusAni, float[] modelEndFuzz){
		return emitVectors(out,window,trueStart,trueStop,varyStart,modelConsensus,table,insideCount,outsideCount,
			model,contigGC,meanLen,offsets,leaveOneOut,stats,writeVectors,featureVersion,locusAni,modelEndFuzz,null,null,null,false);
	}

	static int emitVectors(PrintStream out, byte[] window, int trueStart, int trueStop,
			boolean varyStart, byte[] modelConsensus, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, BaseGraph model, float contigGC, float meanLen,
			int[] offsets, boolean leaveOneOut, VectorStats stats, boolean writeVectors,
			int featureVersion, float locusAni, float[] modelEndFuzz,
			TrnaBoundaryFeatures.NinemerTable familyTable,
			PrintStream trainerOut, String modelName, boolean modelColumn){
		assert(stats!=null) : "G16 wrapper accounting requires a non-null VectorStats accumulator";
		assert(offsets!=null) : "NcrnaFamily boundary offsets define the candidate conservation denominator";
		if(modelColumn&&(trainerOut==null||modelName==null)){throw new IllegalArgumentException("Model-tagged vectors require trainer output and model name");}
		final TrnaBoundaryFeatures.BoundaryType type=(varyStart
			? TrnaBoundaryFeatures.BoundaryType.START : TrnaBoundaryFeatures.BoundaryType.STOP);
		final int trueBoundaryPos=(varyStart ? trueStart : trueStop);
		int n=0;
		for(int offset : offsets){
			int s=(varyStart ? trueStart+offset : trueStart);
			int e=(varyStart ? trueStop : trueStop+offset);
			stats.candidates++;
			if(varyStart){stats.startCandidates++;}else{stats.stopCandidates++;}
			if(offset==0){stats.zeroCandidates++;}
			if(s<0){
				stats.skipStartBeforeWindow++;
				if(offset==0){stats.zeroSkipped++;}
				continue;
			}
			if(e>=window.length){
				stats.skipStopAfterWindow++;
				if(offset==0){stats.zeroSkipped++;}
				continue;
			}
			if(e-s<15){
				stats.skipShort++;
				if(offset==0){stats.zeroSkipped++;}
				continue;
			}
			int label=(offset==0 ? 1 : 0);
			stats.emitted++;
			if(varyStart){stats.startEmitted++;}else{stats.stopEmitted++;}
			if(offset==0){stats.zeroEmitted++;}
			if(!writeVectors){n++; continue;}
			int boundaryPos=(varyStart ? s : e);
			final float ani;
			final float[] fuzz;
			if(featureVersion==NcrnaBoundaryScorer.FEATURES_V2){
				if(!Float.isFinite(locusAni)){throw new IllegalArgumentException("boundaryfeatures=v2 requires finite locus ANI");}
				if(modelEndFuzz==null || modelEndFuzz.length!=3){throw new IllegalArgumentException("boundaryfeatures=v2 requires three model/end fuzz constants");}
				ani=locusAni; fuzz=modelEndFuzz;
			}else if(featureVersion==NcrnaBoundaryScorer.FEATURES_V3){
				if(!Float.isFinite(locusAni)){throw new IllegalArgumentException("boundaryfeatures=v3 requires finite locus identity from QuantumAligner");}
				ani=locusAni;fuzz=null;
			}else{
				final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
				ani=TrnaBoundaryFeatures.aniFeature(candSeq, modelConsensus);
				fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, varyStart);
			}
			float[] prof=(leaveOneOut
				? TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, trueBoundaryPos,
					type, insideCount, outsideCount, table)
				: TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos,
					type, insideCount, outsideCount, table));
			if(featureVersion==NcrnaBoundaryScorer.FEATURES_V3&&familyTable!=null){final float[] familyProf=(leaveOneOut ? TrnaBoundaryFeatures.enrichmentProfile(window,boundaryPos,trueBoundaryPos,type,insideCount,outsideCount,familyTable) : TrnaBoundaryFeatures.enrichmentProfile(window,boundaryPos,type,insideCount,outsideCount,familyTable));blendProfiles(prof,familyProf);}
			float isStop=(varyStart ? 0f : 1f);
			float lengthRatio=(e-s+1)/meanLen;
			if(featureVersion==NcrnaBoundaryScorer.FEATURES_V3){
				if(modelColumn){out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t0.000000\t0.000000\t0.000000\t%s\t%d%n",prof[0],prof[1],prof[2],isStop,lengthRatio,contigGC,ani,modelName,label);}
				else{out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t0.000000\t0.000000\t0.000000\t%d%n",prof[0],prof[1],prof[2],isStop,lengthRatio,contigGC,ani,label);}
				if(trainerOut!=null){trainerOut.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t0.000000\t0.000000\t0.000000\t%d%n",prof[0],prof[1],prof[2],isStop,lengthRatio,contigGC,ani,label);}
			}else{
				if(modelColumn){out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%d%n",ani,prof[0],prof[1],prof[2],isStop,fuzz[0],fuzz[1],fuzz[2],lengthRatio,contigGC,modelName,label);}
				else{out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",ani,prof[0],prof[1],prof[2],isStop,fuzz[0],fuzz[1],fuzz[2],lengthRatio,contigGC,label);}
				if(trainerOut!=null){trainerOut.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",ani,prof[0],prof[1],prof[2],isStop,fuzz[0],fuzz[1],fuzz[2],lengthRatio,contigGC,label);}
			}
			n++;
		}
		return n;
	}

	/** Brian's G28 joint-table rule: equal weight after each table independently
	 * computes its own pseudocount-smoothed profile value. Mutates consensus. */
	static void blendProfiles(float[] consensus,float[] family){
		if(consensus==null||family==null||consensus.length!=family.length){throw new IllegalArgumentException("Family/consensus profile width mismatch");}
		for(int i=0;i<consensus.length;i++){consensus[i]=0.5f*consensus[i]+0.5f*family[i];}
	}

	/** Exact candidate accounting for durable wrapper gates. Whole-locus skips and
	 * per-offset skips are kept separate so a successful record scan cannot be
	 * mistaken for a complete candidate matrix. */
	static final class VectorStats {
		long records, processedLoci;
		long noFlank, malformed, noGC, noModel;
		long candidates, startCandidates, stopCandidates;
		long emitted, startEmitted, stopEmitted;
		long skipStartBeforeWindow, skipStopAfterWindow, skipShort;
		long zeroCandidates, zeroEmitted, zeroSkipped;

		long skipped(){return skipStartBeforeWindow+skipStopAfterWindow+skipShort;}

		void validate(boolean stopOnly, int[] startOffsets, int[] stopOffsets){
			final long wholeSkipped=noFlank+malformed+noGC+noModel;
			if(records!=processedLoci+wholeSkipped){
				throw new IllegalStateException("Record accounting mismatch: records="+records
					+" processed="+processedLoci+" wholeSkipped="+wholeSkipped);
			}
			final int perLocus=(stopOnly ? stopOffsets.length : startOffsets.length+stopOffsets.length);
			final long expectedCandidates=processedLoci*perLocus;
			final long expectedStart=(stopOnly ? 0 : processedLoci*startOffsets.length);
			final long expectedStop=processedLoci*stopOffsets.length;
			if(candidates!=expectedCandidates){
				throw new IllegalStateException("Candidate accounting mismatch: candidates="+candidates
					+" expected="+expectedCandidates+" processed="+processedLoci+" perLocus="+perLocus);
			}
			if(startCandidates!=expectedStart || stopCandidates!=expectedStop){
				throw new IllegalStateException("Per-end candidate accounting mismatch: start="+startCandidates
					+"/"+expectedStart+" stop="+stopCandidates+"/"+expectedStop);
			}
			if(candidates!=emitted+skipped()){
				throw new IllegalStateException("Candidate conservation failed: candidates="+candidates
					+" emitted="+emitted+" skipped="+skipped());
			}
			final long expectedZero=processedLoci*(stopOnly ? 1L : 2L);
			if(zeroCandidates!=expectedZero || zeroEmitted!=expectedZero || zeroSkipped!=0){
				throw new IllegalStateException("Offset-zero gate failed: candidates="+zeroCandidates
					+" emitted="+zeroEmitted+" skipped="+zeroSkipped+" expected="+expectedZero);
			}
		}

		void write(String path, boolean stopOnly, int[] startOffsets, int[] stopOffsets) throws IOException {
			final int perLocus=(stopOnly ? stopOffsets.length : startOffsets.length+stopOffsets.length);
			try(PrintStream out=new PrintStream(new FileOutputStream(path))){
				out.println("metric\tvalue");
				out.println("records\t"+records);
				out.println("processed_loci\t"+processedLoci);
				out.println("whole_locus_skipped\t"+(noFlank+malformed+noGC+noModel));
				out.println("no_flank\t"+noFlank);
				out.println("malformed\t"+malformed);
				out.println("no_gc\t"+noGC);
				out.println("no_model\t"+noModel);
				out.println("candidates_per_locus\t"+perLocus);
				out.println("candidate_total\t"+candidates);
				out.println("start_candidates\t"+startCandidates);
				out.println("stop_candidates\t"+stopCandidates);
				out.println("emitted\t"+emitted);
				out.println("start_emitted\t"+startEmitted);
				out.println("stop_emitted\t"+stopEmitted);
				out.println("skipped_total\t"+skipped());
				out.println("skipped_start_before_window\t"+skipStartBeforeWindow);
				out.println("skipped_stop_after_window\t"+skipStopAfterWindow);
				out.println("skipped_short_candidate\t"+skipShort);
				out.println("zero_candidates\t"+zeroCandidates);
				out.println("zero_emitted\t"+zeroEmitted);
				out.println("zero_skipped\t"+zeroSkipped);
				if(out.checkError()){throw new IOException("Write error: "+path);}
			}
		}
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
