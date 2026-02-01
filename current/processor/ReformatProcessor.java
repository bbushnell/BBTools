package processor;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import shared.Random;
import java.util.concurrent.ConcurrentHashMap;

import cardinality.CardinalityTracker;
import dna.AminoAcid;
import hiseq.IlluminaHeaderParser2;
import hiseq.ReadHeaderParser;
import jgi.CalcTrueQuality;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import stream.FASTQ;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.Quantizer;
import tracker.ReadStats;

/**
 * Encapsulates all read processing logic and statistics for Reformat.
 * Thread-safe via copy() for parallel processing.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 15, 2025
 */
public class ReformatProcessor implements Processor<ReformatProcessor> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public ReformatProcessor(){}

	@Override
	public ReformatProcessor clone(){
		try{
			ReformatProcessor rp=(ReformatProcessor)super.clone();

			//Clone mutable objects that need independent copies
			if(randy!=null){
				rp.randy=shared.Shared.random(sampleseed);
			}
			if(loglog!=null) {
				rp.loglog=loglog.copy();
			}

			return rp;
		}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Parsing              ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean parse(String arg, String a, String b){

		if(a.equals("sample") || a.equals("samplereads") || a.equals("samplereadstarget") || a.equals("srt")){
			sampleReadsTarget=Parse.parseKMG(b);
			sampleReadsExact=(sampleReadsTarget>0);
		}else if(a.equals("samplebases") || a.equals("samplebasestarget") || a.equals("sbt")){
			sampleBasesTarget=Parse.parseKMG(b);
			sampleBasesExact=(sampleBasesTarget>0);
		}else if(a.equals("prioritizelength")){
			prioritizeLength=Parse.parseBoolean(b);
		}else if(a.equals("upsample")){
			allowUpsample=Parse.parseBoolean(b);
		}else if(a.equals("addslash")){
			addslash=Parse.parseBoolean(b);
		}else if(a.equals("addcolon")){
			addcolon=Parse.parseBoolean(b);
		}else if(a.equals("slashspace") || a.equals("spaceslash")){
			boolean x=Parse.parseBoolean(b);
			if(x){
				slash1=" /1";
				slash2=" /2";
			}else{
				slash1="/1";
				slash2="/2";
			}
		}else if(a.equals("addunderscore") || a.equals("underscore")){
			addunderscore=Parse.parseBoolean(b);
		}else if(a.equals("uniquenames")){
			uniqueNames=Parse.parseBoolean(b);
		}else if(a.equals("verifyinterleaved") || a.equals("verifyinterleaving") || a.equals("vint")){
			verifyinterleaving=Parse.parseBoolean(b);
		}else if(a.equals("verifypaired") || a.equals("verifypairing") || a.equals("vpair")){
			verifypairing=Parse.parseBoolean(b);
		}else if(a.equals("allowidenticalnames") || a.equals("ain")){
			allowIdenticalPairNames=Parse.parseBoolean(b);
		}else if(a.equals("rcompmate") || a.equals("rcm")){
			reverseComplementMate=Parse.parseBoolean(b);
		}else if(a.equals("rcomp") || a.equals("rc") || a.equals("reversecomplement")){
			reverseComplement=Parse.parseBoolean(b);
		}else if(a.equals("comp") || a.equals("complement")){
			complement=Parse.parseBoolean(b);
		}else if(a.equals("mappedonly")){
			mappedOnly=Parse.parseBoolean(b);
		}else if(a.equals("pairedonly")){
			pairedOnly=Parse.parseBoolean(b);
		}else if(a.equals("unpairedonly")){
			unpairedOnly=Parse.parseBoolean(b);
		}else if(a.equals("unmappedonly")){
			unmappedOnly=Parse.parseBoolean(b);
		}else if(a.equals("requiredbits") || a.equals("rbits")){
			requiredBits=Parse.parseIntHexDecOctBin(b);
		}else if(a.equals("filterbits") || a.equals("fbits")){
			filterBits=Parse.parseIntHexDecOctBin(b);
		}else if(a.equals("primaryonly")){
			primaryOnly=Parse.parseBoolean(b);
		}else if(a.equals("remap1")){
			remap1=Parse.parseRemap(b);
		}else if(a.equals("remap2")){
			remap2=Parse.parseRemap(b);
		}else if(a.equals("remap")){
			remap1=remap2=Parse.parseRemap(b);
		}else if(a.equals("minmapq")){
			minMapq=Integer.parseInt(b);
		}else if(a.equals("maxmapq")){
			maxMapq=Integer.parseInt(b);
		}else if(a.equals("undefinedton") || a.equals("iupacton") || a.equals("itn")){
			iupacToN=Parse.parseBoolean(b);
		}else if(a.equals("bottom")){
			bottom=Parse.parseBoolean(b);
		}else if(a.equals("top")){
			top=Parse.parseBoolean(b);
		}else if(a.equals("quantize") || a.equals("quantizesticky")){
			quantizeQuality=Quantizer.parse(arg, a, b);
		}else if(a.equals("invert") || a.equals("invertfilters")){
			invertFilters=Parse.parseBoolean(b);
		}else if(a.equals("k")){
			k=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("lastUnderscoreToSpace") || a.equalsIgnoreCase("underscoreToSpace")){
			lastUnderscoreToSpace=Parse.parseBoolean(b);
		}else if(a.equals("fixheader") || a.equals("fixheaders") || a.equals("fixnames")){
			fixHeaders=Parse.parseBoolean(b);
		}else if(a.equals("padleft")){
			padLeft=Integer.parseInt(b);
		}else if(a.equals("padright")){
			padRight=Integer.parseInt(b);
		}else if(a.equals("pad")){
			assert(b!=null) : "No value for pad.";
			if(Character.isLetter(b.charAt(0))){padSymbol=(byte)b.charAt(0);}
			else{padLeft=padRight=Integer.parseInt(b);}
		}else if(a.equals("padsymbol")){
			padSymbol=(byte)b.charAt(0);
		}else if(a.equals("padq")){
			padQ=(byte)Integer.parseInt(b);
		}else if(a.equals("minsubs")){
			minsubs=Integer.parseInt(b);
		}else if(a.equals("tag")){
			tag=b;
		}else if(a.equals("delimiter")){
			delimiter=Parse.parseSymbolToCharacter(b);
			assert(delimiter>0) : delimiter+", '"+b+"'";
		}else if(a.equals("minvalue")){
			minValue=Integer.parseInt(b);
		}else if(a.equals("maxvalue")){
			maxValue=Integer.parseInt(b);
		}else if(a.equals("value")){
			requiredValue=b;
		}else{
			return false;
		}
		return true;
	}

	@Override
	public void setFromParser(Parser parser){
		samplerate=parser.samplerate;
		sampleseed=parser.sampleseed;

		trimBadSequence=parser.trimBadSequence;
		stoptag=SamLine.MAKE_STOP_TAG;

		forceTrimModulo=parser.forceTrimModulo;
		forceTrimLeft=parser.forceTrimLeft;
		forceTrimRight=parser.forceTrimRight;
		forceTrimRight2=parser.forceTrimRight2;
		qtrimLeft=parser.qtrimLeft;
		qtrimRight=parser.qtrimRight;
		trimq=parser.trimq;
		trimE=parser.trimE();
		minAvgQuality=parser.minAvgQuality;
		minAvgQualityBases=parser.minAvgQualityBases;
		chastityFilter=parser.chastityFilter;
		failBadBarcodes=parser.failBadBarcodes;
		removeBadBarcodes=parser.removeBadBarcodes;
		failIfNoBarcode=parser.failIfNoBarcode;
		barcodes=parser.barcodes;
		maxNs=parser.maxNs;
		minConsecutiveBases=parser.minConsecutiveBases;
		minReadLength=parser.minReadLength;
		maxReadLength=parser.maxReadLength;
		minLenFraction=parser.minLenFraction;
		requireBothBad=parser.requireBothBad;
		minGC=parser.minGC;
		maxGC=parser.maxGC;
		filterGC=(minGC>0 || maxGC<1);
		usePairGC=parser.usePairGC;
		tossJunk=parser.tossJunk;
		recalibrateQuality=parser.recalibrateQuality;

		minIdFilter=parser.minIdFilter;
		maxIdFilter=parser.maxIdFilter;
		subfilter=parser.subfilter;
		clipfilter=parser.clipfilter;
		delfilter=parser.delfilter;
		insfilter=parser.insfilter;
		indelfilter=parser.indelfilter;
		dellenfilter=parser.dellenfilter;
		inslenfilter=parser.inslenfilter;
		editfilter=parser.editfilter;

		USE_EDIT_FILTER=(subfilter>-1 || minsubs>-1 || delfilter>-1 || insfilter>-1 || 
			indelfilter>-1 || dellenfilter>-1 || inslenfilter>-1 || editfilter>-1 || clipfilter>-1);

		qtrim=qtrimLeft||qtrimRight;

		if(tag!=null){
			assert(delimiter>0) : 
				"When using a tag, a delimiter must be set; e.g. "
				+ "delimiter=X, delimiter=' ' or delimiter=space\n"
				+ "Most problematic symbols can be spelled out, like tab, "
				+ "pipe, asterisk, greaterthan, etc.";
		}

		if(k>0){parser.loglogk=k;}
		if(parser.loglog && k<1){k=parser.loglogk;}
		loglog=(parser.loglog ? CardinalityTracker.makeTracker(parser) : null);
		pad=padLeft>0 || padRight>0;
	}

	@Override
	public void postParse(){
		pad=(padLeft>0 || padRight>0);
		if(AminoAcid.isFullyDefined(padSymbol)){padQ=Tools.max(padQ, (byte)2);}

		qtrim=qtrimLeft||qtrimRight;

		USE_EDIT_FILTER=(subfilter>-1 || minsubs>-1 || delfilter>-1 || insfilter>-1 || 
			indelfilter>-1 || dellenfilter>-1 || inslenfilter>-1 || editfilter>-1 || clipfilter>-1);

		randy=Shared.threadLocalRandom(sampleseed);

		MAKE_IHIST=ReadStats.COLLECT_INSERT_STATS;
		readstats=ReadStats.collectingStats() ? new ReadStats() : null;

		if(uniqueNames){
			nameMap1=new ObjectIntMap<String>();
			nameMap2=new ObjectIntMap<String>();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Processing Method       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int processReadPair(Read r1, Read r2){

		final int initialLength1=r1.length();
		final int initialLength2=r1.mateLength();

		final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
		final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);

		final SamLine sl1=(r1.samline);

		//Update counters
		readsProcessedT++;
		basesProcessedT+=initialLength1;
		if(r2!=null){
			readsProcessedT++;
			basesProcessedT+=initialLength2;
		}

		//Header fixes
		if(fixHeaders){
			fixHeader(r1);
			fixHeader(r2);
		}

		if(lastUnderscoreToSpace){
			underscoreToSpace(r1);
			underscoreToSpace(r2);
		}

		if(readstats!=null){
			readstats.addToHistograms(r1);

			if(MAKE_IHIST && sl1!=null && !r1.secondary() && sl1.pairnum()==0){
				readstats.addToInsertHistogram(sl1);
			}
		}

		if(loglog!=null){loglog.hash(r1);}

		if(k>0){
			{
				final int kmers=Tools.countKmers(r1.bases, k);
				kmersProcessed+=kmers;
				correctKmers+=(r1.quality==null ? kmers : Tools.countCorrectKmers(r1.quality, k));
			}
			if(r2!=null){
				final int kmers=Tools.countKmers(r2.bases, k);
				kmersProcessed+=kmers;
				correctKmers+=(r2.quality==null ? kmers : Tools.countCorrectKmers(r2.quality, k));
			}
		}

		//Transformations
		if(reverseComplement){r1.reverseComplement();}
		if(complement){r1.complement();}
		if(r2!=null){
			if(reverseComplement || reverseComplementMate){r2.reverseComplementFast();}
			if(complement){r2.complement();}
		}

		//Surface filtering
		if(!bottom || !top){
			rhp.parse(r1);
			if(!bottom && rhp.surface()==1){
				r1.setDiscarded(true);
				if(r2!=null){r2.setDiscarded(true);}
			}
			if(!top && rhp.surface()==2){
				r1.setDiscarded(true);
				if(r2!=null){r2.setDiscarded(true);}
			}
		}

		if(verifypairing){
			String s1=r1==null ? null : r1.id;
			String s2=r2==null ? null : r2.id;
			boolean b=FASTQ.testPairNames(s1, s2, allowIdenticalPairNames);
			if(!b){throw new RuntimeException("Names do not appear to be correctly paired.\n"+s1+"\n"+s2+"\n");}
		}

		//Junk filtering
		if(tossJunk){
			if(r1!=null && r1.junk()){
				lowqBasesT+=r1.length();
				lowqReadsT++;
				r1.setDiscarded(true);
			}
			if(r2!=null && r2.junk()){
				lowqBasesT+=r2.length();
				lowqReadsT++;
				r2.setDiscarded(true);
			}
		}

		//IUPAC to N
		if(iupacToN){
			if(r1!=null){r1.convertUndefinedTo((byte)'N');}
			if(r2!=null){r2.convertUndefinedTo((byte)'N');}
		}

		//Remapping
		if(remap1!=null && r1!=null){
			int swaps=r1.remapAndCount(remap1);
			if(swaps>0){
				basesSwappedT+=swaps;
				readsSwappedT++;
			}
		}
		if(remap2!=null && r2!=null){
			int swaps=r2.remapAndCount(remap2);
			if(swaps>0){
				basesSwappedT+=swaps;
				readsSwappedT++;
			}
		}

		//Trimming
		if(trimBadSequence){
			if(r1!=null){
				int x=TrimRead.trimBadSequence(r1);
				basesQTrimmedT+=x;
				readsQTrimmedT+=(x>0 ? 1 : 0);
			}
			if(r2!=null){
				int x=TrimRead.trimBadSequence(r2);
				basesQTrimmedT+=x;
				readsQTrimmedT+=(x>0 ? 1 : 0);
			}
		}

		//Chastity filter
		if(chastityFilter){
			if(r1!=null && r1.failsChastity()){
				lowqBasesT+=r1.pairLength();
				lowqReadsT+=r1.pairCount();
				r1.setPairDiscarded(true);
			}
		}

		//Barcode filter
		if(removeBadBarcodes){
			if(r1!=null && !r1.discarded() && r1.failsBarcode(barcodes, failIfNoBarcode)){
				lowqBasesT+=r1.pairLength();
				lowqReadsT+=r1.pairCount();
				r1.setPairDiscarded(true);
			}
		}

		//Tag filter
		if(tag!=null && r1!=null && !r1.discarded()){
			boolean pass=passesTagFilter(r1.id, tag, delimiter);
			if(!pass){r1.setPairDiscarded(true);}
		}

		//Bit filters
		if(filterBits!=0 || requiredBits!=0){
			if(r1!=null && !r1.discarded()){
				assert(sl1!=null) : "filterbits and requiredbits only work on sam/bam input.";
				if(((sl1.flag&filterBits)!=0) || ((sl1.flag&requiredBits)!=requiredBits)){
					r1.setDiscarded(true);
					unmappedBasesT+=initialLength1;
					unmappedReadsT++;
				}
			}
		}

		//MAPQ filters
		if(minMapq>=0 || maxMapq>=0){
			if(r1!=null && !r1.discarded()){
				assert(sl1!=null) : "mapq filters only work on sam/bam input.";
				final int mapq=sl1.mapped() ? sl1.mapq : 0;
				if((minMapq>=0 && mapq<minMapq) || (maxMapq>=0 && mapq>maxMapq)){
					r1.setDiscarded(true);
					unmappedBasesT+=initialLength1;
					unmappedReadsT++;
				}
			}
		}

		//ID and edit filters
		if(minIdFilter>=0 || maxIdFilter<1 || USE_EDIT_FILTER){
			if(r1!=null && !r1.discarded()){
				assert(r1.match!=null || r1.samline!=null) : "idfilter requires sam/bam input.";
				boolean pass=passesIDFilter(r1, minIdFilter, maxIdFilter, false);
				if(USE_EDIT_FILTER){
					pass=pass&&passesEditFilter(r1, false);
				}
				if(!pass){
					r1.setDiscarded(true);
					idfilteredBasesT+=initialLength1;
					idfilteredReadsT++;
				}
			}
			if(r2!=null && !r2.discarded()){
				assert(r2.match!=null || r2.samline!=null) : "idfilter requires sam/bam input.";
				boolean pass=passesIDFilter(r2, minIdFilter, maxIdFilter, false);
				if(USE_EDIT_FILTER){
					pass=pass&&passesEditFilter(r2, false);
				}
				if(!pass){
					r2.setDiscarded(true);
					idfilteredBasesT+=initialLength2;
					idfilteredReadsT++;
				}
			}
		}

		if(fixCigar){
			if(SamLine.VERSION==1.3f){
				if(r1!=null && !r1.discarded()){
					assert(sl1!=null) : "Cigar string adjustment only works on sam/bam input.";
					sl1.cigar=SamLine.toCigar13(sl1.cigar);
				}
			}else{
				if(r1!=null && !r1.discarded()){
					assert(sl1!=null) : "Cigar string adjustment only works on sam/bam input.";
					if(r1.match!=null){
						r1.toLongMatchString(false);
						int start=sl1.pos-1;
						int stop=start+Read.calcMatchLength(r1.match)-1;
						sl1.cigar=SamLine.toCigar14(r1.match, start, stop, Integer.MAX_VALUE, r1.bases);
					}
				}
			}
		}

		if(stoptag){
			if(r1!=null && !r1.discarded()){
				assert(sl1!=null) : "stoptag only works on sam/bam input.";
				if(sl1.mapped() && sl1.cigar!=null){
					if(sl1.optional==null){sl1.optional=new ArrayList<String>(2);}
					sl1.optional.add(SamLine.makeStopTag(sl1.pos, sl1.calcCigarLength(false, false), sl1.cigar, r1.perfect()));
				}
			}
		}

		//Paired/unpaired filters
		if(pairedOnly || unpairedOnly){
			assert(sl1!=null) : "pairedonly requires sam/bam input.";
			if(r1!=null && !r1.discarded() && sl1.properPair()!=pairedOnly){
				r1.setDiscarded(true);
				unmappedBasesT+=initialLength1;
				unmappedReadsT++;
			}
		}

		//Mapped/unmapped filters
		if(mappedOnly || unmappedOnly){
			if(r1!=null && !r1.discarded() && (r1.mapped()!=mappedOnly || r1.bases==null || r1.secondary())){
				r1.setDiscarded(true);
				unmappedBasesT+=initialLength1;
				unmappedReadsT++;
			}
			if(r2!=null && !r2.discarded() && (r2.mapped()!=mappedOnly || r2.bases==null || r2.secondary())){
				r2.setDiscarded(true);
				unmappedBasesT+=initialLength2;
				unmappedReadsT++;
			}
		}

		//Primary filter
		if(primaryOnly){
			if(r1!=null && (r1.bases==null || r1.secondary())){
				r1.setDiscarded(true);
				unmappedBasesT+=initialLength1;
				unmappedReadsT++;
			}
			if(r2!=null && (r2.bases==null || r2.secondary())){
				r2.setDiscarded(true);
				unmappedBasesT+=initialLength2;
				unmappedReadsT++;
			}
		}

		//GC filter
		if(filterGC && (initialLength1>0 || initialLength2>0)){
			float gc1=(initialLength1>0 ? r1.gc() : -1);
			float gc2=(initialLength2>0 ? r2.gc() : gc1);
			if(gc1==-1){gc1=gc2;}
			if(usePairGC){
				final float gc;
				if(r2==null){
					gc=gc1;
				}else{
					gc=(gc1*initialLength1+gc2*initialLength2)/(initialLength1+initialLength2);
				}
				gc1=gc2=gc;
			}
			if(r1!=null && !r1.discarded() && (gc1<minGC || gc1>maxGC)){
				r1.setDiscarded(true);
				badGcBasesT+=initialLength1;
				badGcReadsT++;
			}
			if(r2!=null && !r2.discarded() && (gc2<minGC || gc2>maxGC)){
				r2.setDiscarded(true);
				badGcBasesT+=initialLength2;
				badGcReadsT++;
			}
		}

		if(recalibrateQuality){
			if(r1!=null && !r1.discarded()){
				CalcTrueQuality.recalibrate(r1);
			}
			if(r2!=null && !r2.discarded()){
				CalcTrueQuality.recalibrate(r2);
			}
		}

		if(quantizeQuality){
			final byte[] quals1=r1.quality, quals2=(r2==null ? null : r2.quality);
			Quantizer.quantize(quals1);
			Quantizer.quantize(quals2);
		}

		//Force trim
		if(forceTrimLeft>0 || forceTrimRight>=0 || forceTrimModulo>0 || forceTrimRight2>0){
			if(r1!=null && !r1.discarded()){
				final int len=r1.length();
				final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
				final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
				final int b1=forceTrimRight>=0 ? forceTrimRight : len;
				final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
				final int b=Tools.min(b0, b1, b2);
				final int x=TrimRead.trimToPosition(r1, a, b, 1);
				basesFTrimmedT+=x;
				readsFTrimmedT+=(x>0 ? 1 : 0);
				if(r1.length()<minlen1){r1.setDiscarded(true);}
			}
			if(r2!=null && !r2.discarded()){
				final int len=r2.length();
				final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
				final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
				final int b1=forceTrimRight>0 ? forceTrimRight : len;
				final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
				final int b=Tools.min(b0, b1, b2);
				final int x=TrimRead.trimToPosition(r2, a, b, 1);
				basesFTrimmedT+=x;
				readsFTrimmedT+=(x>0 ? 1 : 0);
				if(r2.length()<minlen2){r2.setDiscarded(true);}
			}
		}

		//Quality trim
		if(qtrim){
			if(r1!=null && !r1.discarded()){
				int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1);
				basesQTrimmedT+=x;
				readsQTrimmedT+=(x>0 ? 1 : 0);
			}
			if(r2!=null && !r2.discarded()){
				int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1);
				basesQTrimmedT+=x;
				readsQTrimmedT+=(x>0 ? 1 : 0);
			}
		}

		//Average quality filter
		if(minAvgQuality>0){
			if(r1!=null && !r1.discarded() && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){
				lowqBasesT+=r1.length();
				lowqReadsT++;
				r1.setDiscarded(true);
			}
			if(r2!=null && !r2.discarded() && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
				lowqBasesT+=r2.length();
				lowqReadsT++;
				r2.setDiscarded(true);
			}
		}

		//N filter
		if(maxNs>=0){
			if(r1!=null && !r1.discarded() && r1.countUndefined()>maxNs){
				lowqBasesT+=r1.length();
				lowqReadsT++;
				r1.setDiscarded(true);
			}
			if(r2!=null && !r2.discarded() && r2.countUndefined()>maxNs){
				lowqBasesT+=r2.length();
				lowqReadsT++;
				r2.setDiscarded(true);
			}
		}

		//Consecutive bases filter
		if(minConsecutiveBases>0){
			if(r1!=null && !r1.discarded() && !r1.hasMinConsecutiveBases(minConsecutiveBases)){
				lowqBasesT+=r1.length();
				lowqReadsT++;
				r1.setDiscarded(true);
			}
			if(r2!=null && !r2.discarded() && !r2.hasMinConsecutiveBases(minConsecutiveBases)){
				lowqBasesT+=r2.length();
				lowqReadsT++;
				r2.setDiscarded(true);
			}
		}

		//Length filter
		if(minlen1>0 || minlen2>0 || maxReadLength>0){
			if(r1!=null && !r1.discarded()){
				int rlen=r1.length();
				if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
					r1.setDiscarded(true);
					readShortDiscardsT++;
					baseShortDiscardsT+=rlen;
				}
			}
			if(r2!=null && !r2.discarded()){
				int rlen=r2.length();
				if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
					r2.setDiscarded(true);
					readShortDiscardsT++;
					baseShortDiscardsT+=rlen;
				}
			}
		}

		//Determine if pair should be removed
		boolean remove=false;
		if(r2==null){
			remove=r1.discarded();
		}else{
			remove=requireBothBad ? (r1.discarded() && r2.discarded()) : (r1.discarded() || r2.discarded());
		}
		if(invertFilters){
			remove=!remove;
			r1.setDiscarded(!r1.discarded());
			if(r2!=null){r2.setDiscarded(!r2.discarded());}
		}

		if(remove){return 0;}

		//Padding
		if(pad){
			pad(r1, padLeft, padRight, padSymbol, padQ);
			pad(r2, padLeft, padRight, padSymbol, padQ);
		}

		//Name modifications
		if(uniqueNames || addunderscore || addslash || addcolon){
			if(r1.id==null){r1.id=""+r1.numericID;}
			if(r2!=null && r2.id==null){r2.id=r1.id;}

			if(uniqueNames){
				Integer v=nameMap1.get(r1.id);
				if(v==null){
					nameMap1.put(r1.id, 1);
				}else{
					v++;
					nameMap1.put(r1.id, v);
					r1.id=r1.id+"_"+v;
				}
				if(r2!=null){
					Integer v2=nameMap2.get(r2.id);
					if(v2==null){
						nameMap2.put(r2.id, 1);
					}else{
						v2++;
						nameMap2.put(r2.id, v2);
						r2.id=r2.id+"_"+v2;
					}
				}
			}
			if(addunderscore){
				r1.id=Tools.whitespace.matcher(r1.id).replaceAll("_");
				if(r2!=null){r2.id=Tools.whitespace.matcher(r2.id).replaceAll("_");}
			}
			if(addcolon){
				if(!r1.id.contains(colon1)){r1.id+=colon1;}
				if(r2!=null){
					if(!r2.id.contains(colon2)){r2.id+=colon2;}
				}
			}else if(addslash){
				if(!r1.id.contains(slash1)){r1.id+=slash1;}
				if(r2!=null){
					if(!r2.id.contains(slash2)){r2.id+=slash2;}
				}
			}
		}

		int ret=(!r1.discarded() ? 1 : 0) | (r2!=null && !r2.discarded() ? 2 : 0);
		if(ret==3) {
			pairsOut++;
			pairBasesOut+=r1.pairLength();
		}else if(ret==1){
			singlesOut++;
			singleBasesOut+=r1.length();
		}else if(ret==2){
			singlesOut++;
			singleBasesOut+=r2.length();
		}

		//Fix samline.
		if(!r1.discarded() && r1.samline!=null) {
			SamLine sl=r1.samline;
			sl.seq=r1.bases;
			sl.qual=r1.quality;
			sl.qname=r1.id;
			if(sl.mapped() && sl.strand()==Shared.MINUS) {
				r1.reverseComplementFast();
			}
		}

		return ret;
	}

	@Override
	public boolean processSamLine(SamLine sl){
		Read r=(Read)sl.obj;
		assert(r!=null) : "Input streams need to produce Read objects: "+sl;
		int code=processReadPair(r, null);
		return code>0;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Utility Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int recommendedWorkers() {
		float workers=0.8f;

		// Expensive operations
		if(loglog!=null){workers+=14.5f;}      // Cardinality tracking
		if(recalibrateQuality){workers+=4.0f;}    // Quality recalibration
		if(k>0 && loglog==null){workers+=4.0f;}                 // Kmer counting

		// Moderate operations (0.5f each)
		if(qtrim || qtrimLeft || qtrimRight){workers+=0.5f;}  // Quality trimming
		if(filterGC){workers+=0.5f;}              // GC calculation and filtering
		if(minIdFilter>=0 || maxIdFilter>=0 || USE_EDIT_FILTER){workers+=0.5f;}  // ID/edit filters

		// Light operations (0.25f each)
		if(forceTrimLeft>0 || forceTrimRight>=0 || forceTrimRight2>0 || forceTrimModulo>0){workers+=0.25f;}
		if(minAvgQuality>0 || minAvgQualityBases>0){workers+=0.25f;}
		if(maxNs>=0){workers+=0.25f;}
		if(minConsecutiveBases>0){workers+=0.25f;}
		if(minReadLength>0 || maxReadLength>0){workers+=0.25f;}
		if(reverseComplement || reverseComplementMate || complement){workers+=0.25f;}
		if(trimBadSequence){workers+=0.25f;}
		if(pad){workers+=0.25f;}
		if(remap1!=null || remap2!=null){workers+=0.25f;}
		if(quantizeQuality){workers+=0.25f;}
		if(fixHeaders || lastUnderscoreToSpace || uniqueNames){workers+=0.25f;}
		if(iupacToN){workers+=0.25f;}

		// Very light operations (0.1f each) 
		if(tossJunk || chastityFilter){workers+=0.1f;}
		if(removeBadBarcodes || barcodes!=null){workers+=0.1f;}
		if(mappedOnly || unmappedOnly || pairedOnly || unpairedOnly || primaryOnly){workers+=0.1f;}
		if(filterBits!=0 || requiredBits!=0){workers+=0.1f;}
		if(minMapq>=0 || maxMapq>=0){workers+=0.1f;}
		if(tag!=null){workers+=0.1f;}
		if(minGC>0 || maxGC<1){workers+=0.1f;}
		if(addslash || addcolon || addunderscore){workers+=0.1f;}
		if(verifyinterleaving || verifypairing){workers+=0.1f;}

		// Round up, with minimum 1 and maximum 24
		int result=(int)Math.ceil(workers);
		result=Tools.mid(result, 1, 24);
		if(uniqueNames || sampleReadsExact || sampleBasesExact){result=1;}
		return result;
	}

	public static final void pad(Read r, int padLeft, int padRight, byte padSymbol, byte padQ){
		if(r==null || r.length()==0 || (padLeft<1 && padRight<1)){return;}
		padLeft=Tools.max(0, padLeft);
		padRight=Tools.max(0, padRight);
		r.bases=pad(r.bases, padLeft, padRight, padSymbol);
		r.quality=pad(r.quality, padLeft, padRight, padQ);
	}

	private static final byte[] pad(byte[] old, int padLeft, int padRight, byte padSymbol){
		if(old==null){return null;}
		final int innerLimit=old.length+padLeft;
		byte[] array=new byte[innerLimit+padRight];
		Arrays.fill(array, 0, padLeft, padSymbol);
		Arrays.fill(array, innerLimit, innerLimit+padRight, padSymbol);
		for(int i=0; i<old.length; i++){array[i+padLeft]=old[i];}
		return array;
	}

	public boolean passesTagFilter(String s, String tag, char delimiter){
		if(requiredValue!=null){
			String value=Parse.parseString(s, tag, delimiter);
			return (value==null ? requiredValue==null : value.equals(requiredValue));
		}
		double value=Parse.parseDouble(s, tag, delimiter);
		return value>=minValue && value<=maxValue;
	}

	public static final boolean passesIDFilter(Read r, float minId, float maxId, boolean requireMapped){
		if(!passesMinIDFilter(r, minId, requireMapped)){return false;}
		return passesMaxIDFilter(r, maxId);
	}

	public static final boolean passesMinIDFilter(Read r, float minId, boolean requireMapped){
		if(minId<=0 || r.perfect()){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return !requireMapped;}
		return Read.identityFlat(r.match, true)>=minId;
	}

	public static final boolean passesMaxIDFilter(Read r, float maxId){
		if(maxId>=1){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return true;}
		return Read.identityFlat(r.match, true)<=maxId;
	}

	public final boolean passesEditFilter(Read r, boolean requireMapped){
		if(r.perfect()){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return !requireMapped;}
		r.toLongMatchString(false);

		final int sub=Read.countSubs(r.match);
		final int ins=Read.countInsertions(r.match);
		final int del=Read.countDeletions(r.match);
		final int inscount=Read.countInsertionEvents(r.match);
		final int delcount=Read.countDeletionEvents(r.match);
		final int clip=SamLine.countLeadingClip(r.match)+SamLine.countTrailingClip(r.match);

		boolean bad=false;
		bad=bad||(subfilter>=0 && sub>subfilter);
		bad=bad||(minsubs>=0 && sub<minsubs);
		bad=bad||(clipfilter>=0 && clip>clipfilter);
		bad=bad||(insfilter>=0 && inscount>insfilter);
		bad=bad||(delfilter>=0 && delcount>delfilter);
		bad=bad||(inslenfilter>=0 && r.hasLongInsertion(inslenfilter));
		bad=bad||(dellenfilter>=0 && r.hasLongDeletion(dellenfilter));
		bad=bad||(indelfilter>=0 && inscount+delcount>indelfilter);
		bad=bad||(editfilter>=0 && sub+ins+del>editfilter);

		return !bad;
	}

	public static final void underscoreToSpace(Read r){
		if(r!=null){
			r.id=underscoreToSpace(r.id);
		} 
	}

	public static String underscoreToSpace(String header){
		int x=header.lastIndexOf('_');
		if(x<0){return header;}
		byte[] bytes=header.getBytes();
		bytes[x]=' ';
		return new String(bytes);
	}

	public static final void fixHeader(Read r){
		if(r!=null){
			r.id=fixHeader(r.id);
			if(r.samline!=null){
				r.samline.qname=r.id;
			}
		} 
	}

	public static final String fixHeader(String header){
		if(header==null || header.length()<1){return header;}
		byte[] array=new byte[header.length()];
		boolean changed=false;
		for(int i=0; i<header.length(); i++){
			char c=header.charAt(i);
			byte b=headerSymbols[c];
			array[i]=b;
			changed|=(b!=c);
		}
		if(changed){header=new String(array);}
		return header;
	}

	@Override
	public void add(ReformatProcessor other){
		readsProcessedT+=other.readsProcessedT;
		basesProcessedT+=other.basesProcessedT;

		pairsOut+=other.pairsOut;
		singlesOut+=other.singlesOut;
		pairBasesOut+=other.pairBasesOut;
		singleBasesOut+=other.singleBasesOut;

		basesFTrimmedT+=other.basesFTrimmedT;
		readsFTrimmedT+=other.readsFTrimmedT;

		basesQTrimmedT+=other.basesQTrimmedT;
		readsQTrimmedT+=other.readsQTrimmedT;

		lowqBasesT+=other.lowqBasesT;
		lowqReadsT+=other.lowqReadsT;

		badGcBasesT+=other.badGcBasesT;
		badGcReadsT+=other.badGcReadsT;

		readShortDiscardsT+=other.readShortDiscardsT;
		baseShortDiscardsT+=other.baseShortDiscardsT;

		unmappedReadsT+=other.unmappedReadsT;
		unmappedBasesT+=other.unmappedBasesT;

		idfilteredReadsT+=other.idfilteredReadsT;
		idfilteredBasesT+=other.idfilteredBasesT;

		basesSwappedT+=other.basesSwappedT;
		readsSwappedT+=other.readsSwappedT;

		kmersProcessed+=other.kmersProcessed;
		correctKmers+=other.correctKmers;

		if(loglog!=null) {loglog.add(other.loglog);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void printStats(PrintStream stream) {
		ByteBuilder bb=toStats();
		if(bb.length>0) {bb.nl();}
		stream.print(bb);
	}

	@Override
	public ByteBuilder toStats() {
		ByteBuilder bb=new ByteBuilder();
		final double rpmult=100.0/readsProcessedT, bpmult=100.0/basesProcessedT;

		if(readsSwappedT>0){
			bb.appendln("Base Transforms:        \t"+readsSwappedT+" reads ("+Tools.format("%.2f",readsSwappedT*rpmult)+"%) \t"+
				basesSwappedT+" bases ("+Tools.format("%.2f",basesSwappedT*bpmult)+"%)");
		}
		if(readsQTrimmedT>0 || trimBadSequence){
			bb.appendln("QTrimmed:               \t"+readsQTrimmedT+" reads ("+Tools.format("%.2f",readsQTrimmedT*rpmult)+"%) \t"+
				basesQTrimmedT+" bases ("+Tools.format("%.2f",basesQTrimmedT*bpmult)+"%)");
		}
		if(readsFTrimmedT>0){
			bb.appendln("FTrimmed:               \t"+readsFTrimmedT+" reads ("+Tools.format("%.2f",readsFTrimmedT*rpmult)+"%) \t"+
				basesFTrimmedT+" bases ("+Tools.format("%.2f",basesFTrimmedT*bpmult)+"%)");
		}
		if(readShortDiscardsT>0){
			bb.appendln("Short Read Discards:    \t"+readShortDiscardsT+" reads ("+Tools.format("%.2f",readShortDiscardsT*rpmult)+"%) \t"+
				baseShortDiscardsT+" bases ("+Tools.format("%.2f",baseShortDiscardsT*bpmult)+"%)");
		}
		if(lowqReadsT>0){
			bb.appendln("Low quality discards:   \t"+lowqReadsT+" reads ("+Tools.format("%.2f",lowqReadsT*rpmult)+"%) \t"+
				lowqBasesT+" bases ("+Tools.format("%.2f",lowqBasesT*bpmult)+"%)");
		}
		if(idfilteredReadsT>0){
			bb.appendln("Identity/edit discards: \t"+idfilteredReadsT+" reads ("+Tools.format("%.2f",idfilteredReadsT*rpmult)+"%) \t"+
				idfilteredBasesT+" bases ("+Tools.format("%.2f",idfilteredBasesT*bpmult)+"%)");
		}
		if(badGcReadsT>0){
			bb.appendln("GC content discards:    \t"+badGcReadsT+" reads ("+Tools.format("%.2f",badGcReadsT*rpmult)+"%) \t"+
				badGcBasesT+" bases ("+Tools.format("%.2f",badGcBasesT*bpmult)+"%)");
		}
		if(k>0){
			bb.appendln(k+"-mers processed:        \t"+kmersProcessed);
			bb.appendln("Correct "+k+"-mers:          \t"+Tools.format("%.2f%%", correctKmers*100.0/kmersProcessed));
		}
		if(loglog!=null){
			bb.appendln("Unique "+loglog.k+"-mers:           \t"+loglog.cardinality());
		}
		return bb;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Processing Fields       ----------------*/
	/*--------------------------------------------------------------*/

	//Transformations
	public boolean reverseComplementMate=false;
	public boolean reverseComplement=false;
	public boolean complement=false;
	public boolean iupacToN=false;

	//Verification
	public boolean verifyinterleaving=false;
	public boolean verifypairing=false;
	public boolean allowIdenticalPairNames=true;

	//Filtering
	public boolean tossJunk=false;
	public boolean chastityFilter=false;
	public boolean removeBadBarcodes=false;
	public boolean failBadBarcodes=false;
	public boolean failIfNoBarcode=false;
	public HashSet<String> barcodes=null;

	public boolean mappedOnly=false;
	public boolean pairedOnly=false;
	public boolean unpairedOnly=false;
	public boolean unmappedOnly=false;
	public boolean primaryOnly=false;
	public int requiredBits=0;
	public int filterBits=0;
	public boolean invertFilters=false;

	public boolean bottom=true;
	public boolean top=true;

	//ID/Edit filters
	public float minIdFilter=-1;
	public float maxIdFilter=-1;
	public int subfilter=-1;
	public int minsubs=-1;
	public int clipfilter=-1;
	public int delfilter=-1;
	public int insfilter=-1;
	public int indelfilter=-1;
	public int dellenfilter=-1;
	public int inslenfilter=-1;
	public int editfilter=-1;
	public int minMapq=-1;
	public int maxMapq=-1;
	public boolean USE_EDIT_FILTER=false;

	//Trimming
	public boolean trimBadSequence=false;

	/** Recalibrate quality scores using matrices */
	private boolean recalibrateQuality=false;
	private boolean quantizeQuality=false;
	public boolean qtrim=false;
	public boolean qtrimRight=false;
	public boolean qtrimLeft=false;
	public int forceTrimLeft=0;
	public int forceTrimRight=-1;
	public int forceTrimRight2=0;
	public int forceTrimModulo=0;
	public float trimq=6;
	public float trimE=0.01f;

	//Quality
	public float minAvgQuality=0;
	public int minAvgQualityBases=0;
	public int maxNs=-1;
	public int minConsecutiveBases=0;

	//Length
	public int maxReadLength=0;
	public int minReadLength=0;
	public float minLenFraction=0;
	public boolean requireBothBad=false;

	//GC
	public float minGC=0;
	public float maxGC=1;
	public boolean filterGC=false;
	public boolean usePairGC=false;

	//Header manipulation
	public boolean fixHeaders=false;
	public boolean lastUnderscoreToSpace=false;
	public boolean uniqueNames=false;
	public boolean addslash=false;
	public boolean addcolon=false;
	public boolean addunderscore=false;
	public boolean stoptag=false;


	/** For calculating kmer cardinality */
	public CardinalityTracker loglog;
	public int k=0;

	private ObjectIntMap<String> nameMap1=null;
	private ObjectIntMap<String> nameMap2=null;

	//Padding
	public boolean pad=false;
	public byte padSymbol='N';
	public byte padQ=0;
	public int padLeft=0;
	public int padRight=0;

	//Remapping
	public byte[] remap1=null;
	public byte[] remap2=null;

	//Tag filtering
	public String tag=null;
	public char delimiter=0;
	public float minValue=Float.MIN_VALUE;
	public float maxValue=Float.MAX_VALUE;
	public String requiredValue=null;

	//Sampling
	public float samplerate=1f;
	public long sampleseed=-1;
	public boolean sampleReadsExact=false;
	public boolean sampleBasesExact=false;
	public boolean allowUpsample=false;
	public boolean prioritizeLength=false;
	public long sampleReadsTarget=0;
	public long sampleBasesTarget=0;
	private Random randy=null;

	//Other
	public ReadHeaderParser rhp=new IlluminaHeaderParser2();
	private boolean fixCigar=false;

	/*--------------------------------------------------------------*/
	/*----------------      Statistics Fields       ----------------*/
	/*--------------------------------------------------------------*/

	public boolean MAKE_IHIST;
	public ReadStats readstats;

	public long readsProcessedT=0;
	public long basesProcessedT=0;

	public long pairsOut=0;
	public long singlesOut=0;
	public long pairBasesOut=0;
	public long singleBasesOut=0;

	public long basesFTrimmedT=0;
	public long readsFTrimmedT=0;

	public long basesQTrimmedT=0;
	public long readsQTrimmedT=0;

	public long lowqBasesT=0;
	public long lowqReadsT=0;

	public long badGcBasesT=0;
	public long badGcReadsT=0;

	public long readShortDiscardsT=0;
	public long baseShortDiscardsT=0;

	public long unmappedReadsT=0;
	public long unmappedBasesT=0;

	public long idfilteredReadsT=0;
	public long idfilteredBasesT=0;

	public long basesSwappedT=0;
	public long readsSwappedT=0;

	public long kmersProcessed=0;
	public double correctKmers=0;

	/*--------------------------------------------------------------*/
	/*----------------      Static Fields           ----------------*/
	/*--------------------------------------------------------------*/


	/** For converting headers to filesystem-valid Strings */
	private static final byte[] headerSymbols=new byte[128];

	static{
		Arrays.fill(headerSymbols, (byte)'_');
		for(int i=0; i<128; i++){
			if(Character.isLetterOrDigit(i)){
				headerSymbols[i]=(byte)i;
			}
		}
		char[] acceptable=new char[] {'_', '.', '#', '-', '(', ')', '~'};
		for(char c : acceptable){
			headerSymbols[c]=(byte)c;
		}
	}

	private static String slash1=" /1";
	private static String slash2=" /2";
	private static final String colon1=" 1:";
	private static final String colon2=" 2:";

}