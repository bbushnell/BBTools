package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Streams a roster-complete annotated mitochondrial rRNA census from original GBFF.
 * No product-name selection, SSU classification, pseudogene filtering or sequence deduplication.
 * INSDC Feature Table sections3.4.2.1–2: join concatenates listed segments;
 * complement reverses their order and complements bases. '<'/'>' retain partiality.
 * https://www.insdc.org/submitting-standards/feature-table/
 * Unsupported remote/order/fuzzy locations remain in the ledger, without fabricated sequence.
 * @author Ganyu */
public final class NcrnaMitoRrnaExtract {

	public static void main(String[] args){
		final NcrnaMitoRrnaExtract tool=new NcrnaMitoRrnaExtract(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}
	private NcrnaMitoRrnaExtract(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg:pp.args){
			final String[] parts=Parse.splitOnFirst(arg,'=');
			final String key=parts[0].toLowerCase(java.util.Locale.ROOT),value=parts.length>1?parts[1]:null;
			if(key.equals("roster")){rosterPath=value;}
			else if(key.equals("fna")){fnaPath=value;}
			else if(key.equals("outprefix")){prefix=value;}
			else if(key.equals("expectedrecords")){expected=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		input=parser.in1;
		require(input!=null && rosterPath!=null && prefix!=null && expected>0,"Require in=GBFF roster=TSV outprefix=PATH expectedrecords=N>0");
		require(!parser.append,"Appending would violate the exact roster/copy census");
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
	}
	private void process(){
		loadRoster();
		final String[] suffixes={".candidates.fa.gz",".genomes.fa.gz",".loci.tsv.gz",".records.tsv",".summary.tsv"};
		for(String suffix:suffixes){require(!new File(prefix+suffix).exists(),"Refusing existing census output: "+prefix+suffix);}
		candidates=writer(prefix+suffixes[0]);genomes=writer(prefix+suffixes[1]);loci=writer(prefix+suffixes[2]);records=writer(prefix+suffixes[3]);
		try{
			loci.println("locus_id\taccession\ttaxid\torganelle_class\trrna_index\traw_location\tproduct\tstatus\tpartial\tpseudo\tpseudogene\tstrand\tmin_coordinate\tmax_coordinate\tsegments_in_transcript_order\tcrosses_origin\textracted_length\tgbff");
			records.println("accession\ttaxid\torganelle_class\torganism\tstatus\tbases\ttopology\trrna_features\textracted\tunsupported\tgbff\tsequence_source\tsequence_path\tcontig_location");
			if(fnaPath!=null){readGbff(true);readFallbackFna();}
			readGbff(false);
			int missing=0;
			for(Entry entry:roster){
				if(!entry.seen || !entry.gbffSeen){missing++;}
				b.clear().append(entry.accession).tab().append(entry.taxid).tab();field(b,entry.organelle).tab();field(b,entry.organism).tab();
				b.append(!entry.seen || !entry.gbffSeen?"MISSING":entry.features==0?"NO_RRNA":entry.extracted==0?"UNSUPPORTED_ONLY":"CANDIDATES")
					.tab().append(entry.length).tab().append(entry.topology).tab().append(entry.features).tab().append(entry.extracted).tab().append(entry.features-entry.extracted).tab();
				field(b,input).tab().append(entry.sequenceSource).tab();field(b,entry.sequencePath).tab();
				if(entry.contigLocation.length==0){b.append("NA");}else{field(b,entry.contigLocation);}records.print(b.nl());
			}
			require(missing==0 && seenRecords==expected,"GBFF/roster conservation failed: expected="+expected+" seen="+seenRecords+" missing="+missing);
			require(featureCount==extractedCount+unsupportedCount,"Every annotated rRNA must have exactly one extracted or unsupported ledger row");
			require(fallbackCount>=0 && fallbackCount<=seenRecords,"FASTA fallbacks must be a subset of the complete genome census");
		}finally{close(candidates,genomes,loci,records);}
		final ByteStreamWriter summary=writer(prefix+suffixes[4]);
		try{
			summary.println("roster_records\tgbff_records\tgenome_records\trrna_features\textracted_candidates\tunsupported_locations\tpartial_features\tpseudo_features\tno_rrna_records\ttotal_bases\tfna_fallback_records");
			summary.print(b.clear().append(roster.size()).tab().append(seenRecords).tab().append(seenRecords).tab().append(featureCount)
				.tab().append(extractedCount).tab().append(unsupportedCount).tab().append(partialCount).tab().append(pseudoCount)
				.tab().append(noRrna).tab().append(totalBases).tab().append(fallbackCount).nl());
		}finally{close(summary);}
		outstream.println("MITO_RRNA_EXTRACT_PASS records="+seenRecords+" rrna_features="+featureCount+" extracted="+extractedCount+" unsupported="+unsupportedCount+" fallback_records="+fallbackCount);
	}
	private void loadRoster(){
		final ByteFile bf=ByteFile.makeByteFile(rosterPath,true);final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[] header=bf.nextLine();
			require(header!=null && Tools.equals(header,"accession\ttaxid\torganelle_class\torganism"),"Unexpected mitochondrial roster schema");
			for(byte[] line;(line=bf.nextLine())!=null;){
				lp.set(line);require(lp.terms()==4,"Roster requires four fields");
				final Entry e=new Entry(lp.parseString(0),lp.parseInt(1),lp.parseByteArray(2),lp.parseByteArray(3));
				require(e.accession.matches("[A-Za-z0-9_]+\\.[0-9]+") && e.taxid>0 && e.organelle.length>0,"Require accession.version, positive taxid and organelle class: "+e.accession);
				require(!byAccession.contains(e.accession),"Duplicate roster accession: "+e.accession);
				byAccession.put(e.accession,e);roster.add(e);
			}
		}finally{require(!bf.close(),"Roster input failed: "+rosterPath);}
		require(roster.size()==expected,"Roster count differs from expectedrecords: "+roster.size()+" versus "+expected);
	}
	/** Optional first pass keeps only CONTIG annotation metadata, never all genome sequences. */
	private void readGbff(boolean collectFallback){
		final ByteFile bf=ByteFile.makeByteFile(input,true);Genome g=null;Feature current=null;boolean features=false,origin=false;
		try{
			for(byte[] line;(line=bf.nextLine())!=null;){
				if(Tools.startsWith(line,"LOCUS ")){
					require(g==null,"Missing // terminator before next LOCUS");
					g=new Genome();g.declaredLength=Integer.parseInt(word(line,2));
					require(g.declaredLength>0,"LOCUS must declare a positive length");
					g.circular=hasWord(line,"circular");features=false;origin=false;
				}else if(Tools.startsWith(line,"//")){
					require(g!=null && (origin || g.contig.length()>0),"Record must contain LOCUS and ORIGIN or CONTIG before //");
					if(current!=null){current.finish();}
					finishGenome(g,origin,collectFallback);g=null;current=null;features=false;origin=false;
				}else if(g==null){require(line.length==0,"Unexpected text outside a GBFF record");}
				else if(Tools.startsWith(line,"CONTIG ")){
					if(current!=null){current.finish();}current=null;features=false;origin=false;
					require(g.contig.length()==0,"Duplicate CONTIG header: "+g.accession);
					g.contig.append(trim(line,12));g.contigBlock=true;
				}
				else if(origin){
					for(byte x:line){
						if(x>='a' && x<='z'){x-=32;}
						if(x==' ' || x=='\t' || x>='0' && x<='9'){continue;}
						require("ACGTURYSWKMBDHVN".indexOf(x)>=0,"Invalid ORIGIN base in "+g.accession+": "+(char)x);
						g.bases.append(x);
					}
					require(g.bases.length()<=g.declaredLength,"ORIGIN exceeds LOCUS length for "+g.accession);
				}else if(Tools.startsWith(line,"VERSION ")){
					require(g.accession==null,"Duplicate VERSION line");g.accession=word(line,1);
				}else if(Tools.startsWith(line,"FEATURES ")){features=true;}
				else if(Tools.startsWith(line,"ORIGIN")){
					if(current!=null){current.finish();}current=null;origin=true;features=false;g.contigBlock=false;
				}else if(g.contigBlock && line.length>0 && line[0]==' '){g.contig.append(trim(line,12));
				}else if(features && line.length>5){
					if(line[5]!=' '){
						if(current!=null){current.finish();}current=null;
						if(fieldEquals(line,5,21,"rRNA")){
							current=new Feature();g.features.add(current);current.location.append(trim(line,21));
						}
					}else if(current!=null){current.continuation(line);}
				}
			}
			require(g==null,"Truncated GBFF record at EOF");
		}finally{require(!bf.close(),"GBFF read failed: "+input);}
	}
	private void finishGenome(Genome g,boolean origin,boolean collectFallback){
		require(g.accession!=null,"GBFF record lacks accession.version VERSION");
		final Entry e=byAccession.get(g.accession);require(e!=null,"GBFF accession absent from roster: "+g.accession);
		if(collectFallback){
			require(!e.prescanSeen,"Repeated GBFF accession: "+g.accession);e.prescanSeen=true;e.length=g.declaredLength;
			if(origin){require(g.bases.length()==g.declaredLength,"ORIGIN/LOCUS length mismatch for "+g.accession);}
			else{require(g.bases.length()==0,"Partial ORIGIN cannot be silently replaced: "+g.accession);e.fallback=g;}
			return;
		}
		require(!e.gbffSeen,"Repeated GBFF accession: "+g.accession);e.gbffSeen=true;
		if(origin){emit(g,"GBFF_ORIGIN",input);}
		else{
			require(fnaPath!=null,"CONTIG record requires fna= frozen genomic FASTA: "+g.accession);
			require(e.seen && e.sequenceSource.equals("FNA_CONTIG") && e.length==g.declaredLength && e.features==g.features.size(),"CONTIG record lacks matching completed FNA fallback: "+g.accession);
		}
	}
	private void readFallbackFna(){
		require(fnaPath!=null && new File(fnaPath).isFile(),"Missing frozen fallback FASTA: "+fnaPath);
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(fnaPath,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		int count=0;st.start();
		try{
			for(ListNum<Read> batch;(batch=st.nextList())!=null;){
				for(Read r:batch.list){
					require(r.id!=null,"Fallback FASTA requires accession.version identifiers");
					int end=0;while(end<r.id.length() && !Character.isWhitespace(r.id.charAt(end))){end++;}
					final String accession=r.id.substring(0,end);final Entry e=byAccession.get(accession);
					require(e!=null,"FASTA accession absent from roster: "+accession);
					require(!e.fnaSeen,"Repeated FASTA accession: "+accession);e.fnaSeen=true;count++;
					require(e.prescanSeen,"FASTA accession has no GBFF record: "+accession);
					require(r.bases!=null && r.length()==e.length,"FASTA/LOCUS length mismatch for "+accession);
					if(e.fallback!=null){
						for(byte x:r.bases){require("ACGTURYSWKMBDHVN".indexOf(x)>=0,"Invalid fallback base in "+accession);}
						e.fallback.bases.append(r.bases);emit(e.fallback,"FNA_CONTIG",fnaPath);e.fallback=null;fallbackCount++;
					}
				}st.returnList(batch);
			}
			require(!st.errorState(),"Fallback FASTA stream failed");
		}finally{require(!ReadWrite.closeStream(st),"Fallback FASTA close failed");}
		require(count==expected,"FASTA/roster conservation failed: expected="+expected+" seen="+count);
		for(Entry e:roster){require(e.fnaSeen && e.fallback==null,"Missing fallback FASTA record: "+e.accession);}
	}
	private void emit(Genome g,String sequenceSource,String sequencePath){
		assert(g!=null):"Only complete GBFF records enter emission";
		require(g.accession!=null,"GBFF record lacks accession.version VERSION");
		final Entry e=byAccession.get(g.accession);
		require(e!=null,"GBFF accession absent from roster: "+g.accession);
		require(!e.seen,"Repeated GBFF accession: "+g.accession);
		require(g.bases.length()==g.declaredLength,"ORIGIN/LOCUS length mismatch for "+g.accession+": "+g.bases.length()+" vs "+g.declaredLength);
		e.seen=true;e.length=g.declaredLength;e.topology=g.circular?"circular":"linear_or_unspecified";e.features=g.features.size();
		e.sequenceSource=sequenceSource;e.sequencePath=sequencePath;e.contigLocation=g.contig.toBytes();
		seenRecords++;totalBases+=e.length;if(e.features==0){noRrna++;}
		b.clear().append(e.accession).append(" tid_").append(e.taxid);fasta(genomes,b.toString(),g.bases.array,g.bases.length());
		int ordinal=0;
		for(Feature f:g.features){
			ordinal++;featureCount++;
			final String raw=new String(f.location.array,0,f.location.length(),StandardCharsets.UTF_8),id=e.accession+"_rrna_"+ordinal;
			final boolean partial=f.partial || raw.indexOf('<')>=0 || raw.indexOf('>')>=0;
			if(partial){partialCount++;}if(f.pseudo){pseudoCount++;}
			String status="EXTRACTED",strand="NA",segments="NA";int lo=0,hi=0;boolean crosses=false;
			sequence.clear();
			try{
				// TODO: GbffFeature.parseStartStop collapses joins, loses partial markers and
				// rejects origin-spanning start>stop (constructor), then GbffLocus drops errors.
				// Read original locations here; never use that derived GFF as the full census.
				final ArrayList<Segment> spans=new Location(raw,e.length).parse();
				final ByteBuilder geometry=new ByteBuilder();int direction=0,previous=0;lo=Integer.MAX_VALUE;
				for(Segment s:spans){
					lo=Math.min(lo,s.lo);hi=Math.max(hi,s.hi);
					if(direction==0){direction=s.strand;}else if(direction!=s.strand){direction=2;}
					final int first=s.strand==1?s.lo:s.hi;
					if(previous>0 && g.circular && (s.strand==1 && first<previous || s.strand==-1 && first>previous)){crosses=true;}
					previous=s.strand==1?s.hi:s.lo;
					if(geometry.length()>0){geometry.append(',');}geometry.append(s.lo).append("..").append(s.hi).append(s.strand==1?":+":":-");
					if(s.strand==1){sequence.append(g.bases.array,s.lo-1,s.hi-s.lo+1);}
					else{for(int i=s.hi-1;i>=s.lo-1;i--){sequence.append(AminoAcid.baseToComplementExtended[g.bases.array[i]]);}}
				}
				strand=direction==1?"+":direction==-1?"-":"MIXED";segments=geometry.toString();
				require(sequence.length()>0,"A supported rRNA location must extract bases: "+id);
				b.clear().append(id).append(" tid_").append(e.taxid).append(" source_accession=").append(e.accession)
					.append(" partial=").append(partial).append(" pseudo=").append(f.pseudo);
				fasta(candidates,b.toString(),sequence.array,sequence.length());e.extracted++;extractedCount++;
			}catch(UnsupportedLocation ex){status=ex.getMessage();unsupportedCount++;}
			b.clear().append(id).tab().append(e.accession).tab().append(e.taxid).tab();field(b,e.organelle).tab().append(ordinal).tab();
			field(b,f.location).tab();field(b,f.product).tab().append(status).tab().append(partial).tab().append(f.pseudo).tab();
			field(b,f.pseudogene).tab().append(strand).tab().append(lo).tab().append(hi).tab()
				.append(segments).tab().append(crosses).tab().append(sequence.length()).tab();field(b,input).nl();loci.print(b);
		}
	}
	/** Limited, explicit INSDC local-location grammar. Syntax errors in supported forms fail loudly. */
	private static final class Location {
		Location(String text_,int length_){text=text_;length=length_;}
		ArrayList<Segment> parse(){
			require(!text.isEmpty(),"Empty rRNA location");
			if(text.indexOf(':')>=0){throw new UnsupportedLocation("UNSUPPORTED_REMOTE");}
			if(text.contains("order(")){throw new UnsupportedLocation("UNSUPPORTED_ORDER");}
			if(text.indexOf('^')>=0 || text.indexOf('?')>=0 || text.contains("one-of") || text.contains("within")){throw new UnsupportedLocation("UNSUPPORTED_FUZZY");}
			for(int i=0;i<text.length();i++){
				if(text.charAt(i)=='.' && (i==0 || text.charAt(i-1)!='.') && (i+1==text.length() || text.charAt(i+1)!='.')){
					throw new UnsupportedLocation("UNSUPPORTED_FUZZY");
				}
			}
			expr(0);space();require(pos==text.length(),"Trailing location syntax: "+text);
			require(!spans.isEmpty(),"Location contains no spans: "+text);return spans;
		}
		void expr(int depth){
			require(depth<64,"Location nesting exceeds64: "+text);space();
			if(take("complement(")){
				final int first=spans.size();expr(depth+1);expect(')');
				Collections.reverse(spans.subList(first,spans.size()));
				for(int i=first;i<spans.size();i++){spans.get(i).strand=-spans.get(i).strand;}
			}else if(take("join(")){
				expr(depth+1);int members=1;space();
				while(pos<text.length() && text.charAt(pos)==','){pos++;expr(depth+1);members++;space();}
				require(members>=2,"join requires multiple segments: "+text);expect(')');
			}else{
				if(pos<text.length() && (Character.isLetter(text.charAt(pos)) || text.charAt(pos)=='(')){throw new UnsupportedLocation("UNSUPPORTED_OPERATOR_OR_FUZZY");}
				final int lo=coordinate();space();int hi=lo;
				if(take("..")){hi=coordinate();}
				else if(pos<text.length() && text.charAt(pos)=='.'){throw new UnsupportedLocation("UNSUPPORTED_FUZZY");}
				require(lo>=1 && hi>=lo && hi<=length,"Location outside source sequence or reversed interval: "+text+" length="+length);
				spans.add(new Segment(lo,hi));
			}
		}
		int coordinate(){
			space();if(pos<text.length() && (text.charAt(pos)=='<' || text.charAt(pos)=='>')){pos++;}
			if(pos<text.length() && text.charAt(pos)=='('){throw new UnsupportedLocation("UNSUPPORTED_FUZZY");}
			final int first=pos;long n=0;
			while(pos<text.length() && text.charAt(pos)>='0' && text.charAt(pos)<='9'){
				n=n*10+text.charAt(pos++)-'0';require(n<=Integer.MAX_VALUE,"Location coordinate overflow: "+text);
			}
			require(pos>first,"Missing coordinate: "+text);return (int)n;
		}
		boolean take(String token){assert(token!=null && !token.isEmpty()):"Grammar tokens must advance the location parser";space();if(!text.startsWith(token,pos)){return false;}pos+=token.length();return true;}
		void expect(char c){space();require(pos<text.length() && text.charAt(pos)==c,"Expected "+c+" in location: "+text);pos++;}
		void space(){assert(pos>=0 && pos<=text.length()):"Location cursor cannot move outside its original text";while(pos<text.length() && Character.isWhitespace(text.charAt(pos))){pos++;}}
		final String text;final int length;int pos;final ArrayList<Segment> spans=new ArrayList<Segment>();
	}
	private static final class Segment{Segment(int a,int z){lo=a;hi=z;}final int lo,hi;int strand=1;}
	private static final class UnsupportedLocation extends RuntimeException{UnsupportedLocation(String reason){super(reason);}private static final long serialVersionUID=1L;}
	private static final class Feature {
		void continuation(byte[] line){
			assert(line!=null):"ByteFile EOF is handled before feature continuation";
			final byte[] value=trim(line,21);
			if(value.length>0 && value[0]=='/'){flushQualifier();qualifier.append(value);qualifiers=true;}
			else if(qualifiers){if(qualifier.length()>0){qualifier.append(' ');}qualifier.append(value);}
			else{location.append(value);}
		}
		void finish(){flushQualifier();require(location.length()>0,"rRNA feature missing location");}
		void flushQualifier(){
			assert(qualifier.length()==0 || qualifier.array[0]=='/'):"Qualifier continuation must retain its original slash-prefixed key";
			if(qualifier.startsWith("/product=")){if(product.length()>0){product.append(" | ");}product.append(qualifier.array,9,qualifier.length()-9);}
			else if(qualifier.startsWith("/pseudogene=")){pseudo=true;pseudogene.append(qualifier.array,12,qualifier.length()-12);}
			else if(fieldEquals(qualifier.array,0,qualifier.length(),"/pseudo")){pseudo=true;}
			else if(fieldEquals(qualifier.array,0,qualifier.length(),"/partial")){partial=true;}
			qualifier.clear();
		}
		final ByteBuilder location=new ByteBuilder(),product=new ByteBuilder(),pseudogene=new ByteBuilder(),qualifier=new ByteBuilder();
		boolean qualifiers,pseudo,partial;
	}
	private static final class Genome {String accession;int declaredLength;boolean circular,contigBlock;final ByteBuilder bases=new ByteBuilder(),contig=new ByteBuilder();final ArrayList<Feature> features=new ArrayList<Feature>();}
	private static final class Entry {
		Entry(String a,int t,byte[] c,byte[] o){accession=a;taxid=t;organelle=c;organism=o;}
		final String accession;final byte[] organelle,organism;final int taxid;boolean seen,gbffSeen,prescanSeen,fnaSeen;
		int length,features,extracted;Genome fallback;byte[] contigLocation=new byte[0];String topology="NA",sequenceSource="NA",sequencePath="NA";
	}
	private static String word(byte[] line,int index){
		int n=0;for(int i=0;i<line.length;){while(i<line.length && asciiSpace(line[i])){i++;}final int from=i;while(i<line.length && !asciiSpace(line[i])){i++;}
			if(from<i && n++==index){return new String(line,from,i-from,StandardCharsets.US_ASCII);}}
		throw new IllegalArgumentException("Missing field "+index+" in GBFF header");
	}
	private static boolean hasWord(byte[] line,String value){
		assert(line!=null && value!=null):"LOCUS token checks receive an actual input line and fixed token";
		for(int i=0;i<line.length;){while(i<line.length && asciiSpace(line[i])){i++;}final int from=i;while(i<line.length && !asciiSpace(line[i])){i++;}
			if(fieldEquals(line,from,i,value)){return true;}}return false;
	}
	private static boolean fieldEquals(byte[] line,int from,int to,String value){
		assert(line!=null && from>=0 && to>=from && value!=null):"Feature/header field ranges must be ordered before clipping";
		to=Math.min(to,line.length);while(from<to && asciiSpace(line[from])){from++;}while(to>from && asciiSpace(line[to-1])){to--;}
		if(to-from!=value.length()){return false;}for(int i=0;i<value.length();i++){if(line[from+i]!=value.charAt(i)){return false;}}return true;
	}
	private static byte[] trim(byte[] line,int from){
		assert(line!=null && from>=0):"GBFF field starts are nonnegative fixed columns";
		int to=line.length;from=Math.min(from,to);while(from<to && asciiSpace(line[from])){from++;}while(to>from && asciiSpace(line[to-1])){to--;}
		return java.util.Arrays.copyOfRange(line,from,to);
	}
	private static ByteBuilder field(ByteBuilder bb,String text){
		assert(text!=null):"Census fields use explicit NA rather than null";
		return field(bb,text.getBytes(StandardCharsets.UTF_8));
	}
	private static ByteBuilder field(ByteBuilder bb,ByteBuilder text){return text.length()==0?bb.append("NA"):field(bb,text.array,text.length());}
	private static ByteBuilder field(ByteBuilder bb,byte[] text){return field(bb,text,text.length);}
	/** Escape only TSV controls and percent; every other original byte, including UTF-8, is unchanged. */
	private static ByteBuilder field(ByteBuilder bb,byte[] text,int length){
		assert(length>=0 && length<=text.length):"Retained metadata bounds must fit its original byte buffer";
		for(int i=0;i<length;i++){final byte c=text[i];if(c=='%' || c=='\t' || c=='\r' || c=='\n'){
			bb.append('%').append("0123456789ABCDEF".charAt(c>>4)).append("0123456789ABCDEF".charAt(c&15));
		}else{bb.append(c);}}return bb;
	}
	private static boolean asciiSpace(byte c){return c>=0 && c<=' ';}
	private void fasta(ByteStreamWriter out,String header,byte[] bases,int length){
		assert(length>0 && length<=bases.length):"FASTA records contain validated nonempty sequence";
		out.print(b.clear().append('>').append(header).nl());
		for(int i=0;i<length;i+=1023){out.print(b.clear().append(bases,i,Math.min(1023,length-i)).nl());}
	}
	private static ByteStreamWriter writer(String path){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;}
	private static void close(ByteStreamWriter... streams){boolean failed=false;for(ByteStreamWriter w:streams){failed|=w.poisonAndWait();}require(!failed,"Census output writer failed");}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private String input,rosterPath,prefix,fnaPath;private int expected,seenRecords,noRrna,fallbackCount;
	private long featureCount,extractedCount,unsupportedCount,partialCount,pseudoCount,totalBases;
	private final PrintStream outstream;
	private final ObjectMap<String,Entry> byAccession=new ObjectMap<String,Entry>(String.class,Entry.class);
	private final ArrayList<Entry> roster=new ArrayList<Entry>();
	private final ByteBuilder b=new ByteBuilder(),sequence=new ByteBuilder();
	private ByteStreamWriter candidates,genomes,loci,records;
}
