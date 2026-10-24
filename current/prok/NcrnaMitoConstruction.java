package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;

/** Audits accepted mitochondrial CM evidence and exports the declared SSU-only
 * construction cohorts without filtering partial, pseudo or discordant members.
 * Raw-bit-score ranking and short-length categories are descriptive, not ancestry
 * or a second admission rule. CM-only exact-coordinate groups are audit-only.
 * @author Brian, Raiden */
public final class NcrnaMitoConstruction {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoConstruction.class,false);final Parser parser=new Parser();
		final NcrnaMitoConstruction x=new NcrnaMitoConstruction();String verified=null;
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("loci")){x.locusPath=v;}else if(k.equals("hits")){x.hitPath=v;}else if(k.equals("associations")){x.associationPath=v;}
				else if(k.equals("pairs")){x.pairPath=v;}else if(k.equals("genomes")){x.genomePath=v;}
				else if(k.equals("verified")){verified=v;}else if(k.equals("outprefix")){x.prefix=v;}
				else if(k.equals("expectedloci")){x.expectedLoci=Parse.parseIntKMG(v);}else if(k.equals("expectedgenomes")){x.expectedGenomes=Parse.parseIntKMG(v);}
				else if(k.equals("expectedhits")){x.expectedHits=Parse.parseIntKMG(v);}else if(k.equals("expectedgroups")){x.expectedGroups=Parse.parseIntKMG(v);}
				else if(k.equals("expected16")){x.expected16=Parse.parseIntKMG(v);}else if(k.equals("expected18")){x.expected18=Parse.parseIntKMG(v);}
				else if(k.equals("expectedunsupported")){x.expectedUnsupported=Parse.parseIntKMG(v);}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			x.input=parser.in1;
			require(verified!=null && x.prefix!=null && !parser.append && parser.maxReads<0,"Require verified inputs, fresh outprefix, and the complete census without append/read limits");
			require(x.expectedLoci>0 && x.expectedGenomes>0 && x.expectedHits>0 && x.expectedGroups>0 && x.expected16>=0 && x.expected18>=0 && x.expectedUnsupported>=0,"Positive source counts and explicit nonnegative cohort counts protect the declared full population");
			final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);
			for(String path:new String[]{x.input,x.locusPath,x.hitPath,x.associationPath,x.pairPath,x.genomePath}){require(path!=null && pins.get(path)!=null,"Every construction input must be bound by the producer: "+path);}
			Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
			RrnaResourceIO.fresh(x.prefix,SUFFIXES);x.loadGenomes();x.loadLoci();x.loadPairs();x.loadHits();x.loadAssociations();x.classify();x.loadSequences();x.write();
			pp.outstream.println("MITO_CONSTRUCTION_PASS loci="+x.loci.size()+" m16S="+x.cohortCounts[0]+" m18S="+x.cohortCounts[1]+" unsupported="+x.cohortCounts[4]+" cm_only_exact_groups="+x.cmOnly.size()+" short_N="+x.shortN+" no_production_change=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void loadGenomes(){
		read(genomePath,9,"accession",p->{final String id=p.parseString(0);final int tid=p.parseInt(1);
			require(tid>0 && taxids.get(id)<0,"Genome keys must be unique and retain positive host taxids");taxids.put(id,tid);
		});require(taxids.size()==expectedGenomes,"Audit genome population must equal the accepted evidence census");
	}
	private void loadLoci(){
		read(locusPath,28,"locus_id",p->{
			final Locus l=new Locus();l.id=p.parseString(0);l.accession=p.parseString(1);l.uniqueId=p.parseString(2);l.status=p.parseString(8);
			l.taxid=p.parseInt(11);l.product=p.parseString(15);l.partial=bit(p,17);l.pseudo=bit(p,18);l.segments=p.parseString(23);l.crosses=bit(p,24);l.length=p.parseInt(25);
			require(p.termEquals(l.id,9) && p.termEquals(l.accession,10) && l.taxid==taxids.get(l.accession) && p.termEquals("EXTRACTED",16) && l.length>0,"Audit rows must preserve exact extracted source identity, host and length");
			require(byLocus.put(l.id,l)==null,"Duplicate annotated locus: "+l.id);loci.add(l);
			Unique u=byUnique.get(l.uniqueId);if(u==null){u=new Unique();byUnique.put(l.uniqueId,u);}u.copies++;l.unique=u;
		});require(loci.size()==expectedLoci,"Every accepted source annotation must remain in the audit");
	}
	private void loadPairs(){
		assert(loci.size()==expectedLoci):"Universal comparisons annotate the complete source census loaded first";
		final int[] n={0};read(pairPath,18,"locus_id",p->{final Locus l=byLocus.get(p.parseString(0));
			require(l!=null && p.termEquals(l.uniqueId,1) && p.parseInt(3)==l.length,"Universal comparison must bind the same original sequence");
			final int side=p.termEquals("universal",2)?0:p.termEquals("euk18S_universal",2)?1:-1;
			require(side>=0 && Double.isNaN(l.identity[side]),"Each locus requires exactly one comparison per universal model");
			final double id=p.parseDouble(6);require(Double.isFinite(id) && id>=0 && id<=1,"Universal identity must be a finite fraction");l.identity[side]=id;n[0]++;
		});require(n[0]==2*expectedLoci,"Both universal comparisons must cover every annotation");
	}
	private void loadHits(){
		hits=new Hit[expectedGroups];final int[] count={0,0};
		read(hitPath,28,"run_id",p->{count[0]++;final int group=p.parseInt(20);require(group>=0 && group<hits.length,"Physical group index must fit the accepted complete roster");
			if(!bit(p,21)){return;}final Hit h=new Hit();h.group=group;h.run=p.parseString(0);h.line=p.parseLong(2);h.model=p.parseString(3);h.ssu=modelIndex(h.model);
			h.genomic=p.termEquals("genomic",5);require(h.genomic || p.termEquals("candidate",5),"Unknown target universe");
			h.accession=p.parseString(7);h.physical=p.parseString(17);h.crosses=bit(p,18);h.trunc=p.parseString(22);h.bits=p.parseDouble(24);h.coverage=p.parseDouble(12);
			h.length=physicalLength(h.physical);
			require(h.line>0 && taxids.get(h.accession)>0 && Double.isFinite(h.bits) && h.coverage>0 && h.coverage<=1 && h.length>0,"Representative hits need valid source coordinates and model coverage");
			require((h.ssu>=0)==p.termEquals("SSU",4) && hits[group]==null,"Each exact group has one representative with its declared family");hits[group]=h;count[1]++;
		});require(count[0]==expectedHits && count[1]==expectedGroups,"Raw and representative counts must equal the accepted reduction");
	}
	private void loadAssociations(){
		assert(hits!=null && hits.length==expectedGroups):"Associations resolve against the complete accepted representative roster";
		read(associationPath,9,"run_id",p->{if(!bit(p,3)){return;}final int g=p.parseInt(2);require(g>=0 && g<hits.length,"Association names an undeclared physical group");
			final Hit h=hits[g];final Locus l=byLocus.get(p.parseString(4));final int same=p.parseInt(5),opposite=p.parseInt(6);
			require(h!=null && l!=null && l.accession.equals(h.accession) && p.termEquals(h.run,0) && p.parseLong(1)==h.line,"Association must name the exact representative and source accession");
			require(p.parseInt(7)==l.length && p.parseInt(8)==h.length && same>=0 && opposite>=0 && same+opposite>0 && same+opposite<=Math.min(l.length,h.length),"Overlap accounting must conserve the actual source/hit footprints");
			l.anyHit=true;if(same>0){h.annotated=true;if(h.ssu>=0){l.best.offer(h,same);}else{l.lsu=true;}}
		});
		final ObjectMap<String,CmOnly> map=new ObjectMap<String,CmOnly>(String.class,CmOnly.class);
		for(Hit h:hits){require(h!=null,"Missing representative cannot become absent evidence");if(h.genomic && h.ssu>=0 && !h.annotated){
			final String key=h.accession+'\t'+h.physical;CmOnly c=map.get(key);if(c==null){c=new CmOnly();c.accession=h.accession;c.physical=h.physical;map.put(key,c);cmOnly.add(c);}c.members++;c.best.offer(h,0);
		}}
	}
	private void classify(){
		shortN=Integer.MAX_VALUE;
		for(Locus l:loci){final boolean ssu=l.best.winner()>=0;
			final String observed=!ssu && !l.lsu?"NO_SAME_STRAND_GA_HIT_UNRESOLVED":ssu && l.lsu?"SSU_LSU_CONFLICT_REVIEW":ssu?"SSU_SUPPORT_REVIEW":"LSU_SUPPORT_REVIEW";
			require(observed.equals(l.status),"Recomputed evidence support must agree with frozen classification: "+l.id);
			require(!Double.isNaN(l.identity[0]) && !Double.isNaN(l.identity[1]),"Missing paired universal identity");
			l.preference=l.identity[0]>l.identity[1]?0:l.identity[0]==l.identity[1]?1:2;
			l.cohort=ssu && !l.lsu?(l.preference==0?0:l.preference==2?1:2):ssu?3:!l.lsu?4:5;
			cohortCounts[l.cohort]++;l.unique.cohortCopies[l.cohort]++;
			if(ssu && !l.partial && !l.pseudo){shortN=Math.min(shortN,l.length);}
		}
		require(shortN<Integer.MAX_VALUE,"Descriptive N requires at least one non-partial/non-pseudo annotation with same-strand SSU support");
		require(cohortCounts[0]==expected16 && cohortCounts[1]==expected18 && cohortCounts[4]==expectedUnsupported,"Frozen construction cohorts cannot change through audit flags or CM ranking");
		for(Locus l:loci){if(l.cohort==4){unsupported[unsupportedCategory(l.product,l.length,shortN,l.anyHit)]++;}}
	}
	private void loadSequences(){
		assert(shortN>0 && cohortCounts[0]==expected16 && cohortCounts[1]==expected18):"Sequence export follows the locked evidence/preference classification, not a new sequence filter";
		final int[] n={0};RrnaResourceIO.stream(input,r->{final Locus l=byLocus.get(RrnaResourceIO.first(r.id));
			require(l!=null && !l.sequenceSeen && r.bases!=null && r.length()==l.length,"FASTA must contain each exact original annotation once, without length changes");
			l.sequenceSeen=true;n[0]++;if(l.unique.bases==null){l.unique.bases=r.bases.clone();}else{require(Arrays.equals(l.unique.bases,r.bases),"A shared unique-sequence ID must bind identical actual bases");}
		});require(n[0]==expectedLoci,"Export input must conserve the complete annotated corpus, not a benchmark subset");
		for(Locus l:loci){require(l.sequenceSeen,"Missing source FASTA member: "+l.id);}
	}
	private void write(){
		assert(loci.size()==expectedLoci):"Audit output must cover the original full annotation denominator";
		final ByteStreamWriter audit=RrnaResourceIO.writer(prefix+SUFFIXES[0]),cross=RrnaResourceIO.writer(prefix+SUFFIXES[1]),reasons=RrnaResourceIO.writer(prefix+SUFFIXES[2]),summary=RrnaResourceIO.writer(prefix+SUFFIXES[3]);
		final ByteBuilder b=new ByteBuilder();final int[][][] counts=new int[3][5][5];
		try{
			audit.println("kind\tid\taccession\ttaxid\tconstruction\tcm_status\tlength\tpartial\tpseudo\tproduct\tsource_segments\tsource_crosses_origin\tunique_sequence_id\tsource_copies\tconstruction_copies\tuniversal_preference\tquantum_16S\tquantum_18S\tbest_ssu_cm\tbest_ssu_bits\tmodel_coverage\ttruncation\tbest_tied_models\tcm_concordance\tsame_strand_overlap_bp\tsource_overlap_fraction\thit_overlap_fraction\thit_segments\texact_genomic_evidence_members\tunsupported_description");
			for(Locus l:loci){b.clear().append("ANNOTATION").tab().append(l.id).tab().append(l.accession).tab().append(l.taxid).tab().append(COHORTS[l.cohort]).tab().append(l.status).tab().append(l.length).tab().append(l.partial).tab().append(l.pseudo).tab().append(l.product).tab().append(l.segments).tab().append(l.crosses).tab().append(l.uniqueId).tab().append(l.unique.copies).tab().append(l.unique.cohortCopies[l.cohort]).tab().append(PREFERENCES[l.preference]).tab().append(l.identity[0],9).tab().append(l.identity[1],9);
				appendChoice(b,l.best,l.preference,l.length);b.append("\tNA\t").append(l.cohort==4?REASONS[unsupportedCategory(l.product,l.length,shortN,l.anyHit)]:"NA");audit.print(b.nl());
				if(l.cohort<=2){counts[l.preference][l.best.winner()][l.best.concordance(l.preference)]++;}
			}
			int index=0;for(CmOnly c:cmOnly){final Hit h=c.best.best[c.best.winner()];
				b.clear().append("CM_ONLY_EXACT_GEOMETRY\tCM_ONLY_").append(index++).tab().append(c.accession).tab().append(taxids.get(c.accession)).append("\tAUDIT_ONLY\tSSU_CM_ONLY_REVIEW\t").append(h.length).append("\tNA\tNA\tNA\t").append(c.physical).tab().append(h.crosses).append("\tNA\tNA\tNA\tNA\tNA\tNA");
				appendChoice(b,c.best,-1,-1);audit.print(b.tab().append(c.members).append("\tNA").nl());
			}
			cross.println("universal_preference\thighest_raw_bit_score_ssu_cm\tannotations\tconcordant\tdiscordant\tspecialized_mito_model\ttied_best_models\tuniversal_tie");
			for(int p=0;p<3;p++){for(int m=0;m<5;m++){int n=0;for(int c:counts[p][m]){n+=c;}b.clear().append(PREFERENCES[p]).tab().append(MODELS[m]).tab().append(n);for(int c:counts[p][m]){b.tab().append(c);}cross.print(b.nl());}}
			reasons.println("precedence\tdescriptive_category\tannotations\tlength_N\tselection_effect");
			for(int i=0;i<REASONS.length;i++){reasons.print(b.clear().append(i+1).tab().append(REASONS[i]).tab().append(unsupported[i]).tab().append(shortN).append("\tNONE").nl());}
			summary.println("construction\tcopies\tunique_sequences\tnonpartial_nonpseudo\tpartial\tpseudo\tmin_length\tmax_length\tlength_N");
			for(int group=0;group<COHORTS.length;group++){int unique=0,intact=0,partial=0,pseudo=0,min=Integer.MAX_VALUE,max=0;final ObjectIntMap<String> seen=new ObjectIntMap<String>(String.class);
				for(Locus l:loci){if(l.cohort!=group){continue;}if(seen.get(l.uniqueId)<0){seen.put(l.uniqueId,1);unique++;}if(!l.partial && !l.pseudo){intact++;}if(l.partial){partial++;}if(l.pseudo){pseudo++;}min=Math.min(min,l.length);max=Math.max(max,l.length);}
				summary.print(b.clear().append(COHORTS[group]).tab().append(cohortCounts[group]).tab().append(unique).tab().append(intact).tab().append(partial).tab().append(pseudo).tab().append(min==Integer.MAX_VALUE?0:min).tab().append(max).tab().append(shortN).nl());
			}
		}finally{RrnaResourceIO.finish(audit,cross,reasons,summary);}
		writeFasta(0,SUFFIXES[4]);writeFasta(1,SUFFIXES[5]);
	}
	private void writeFasta(int cohort,String suffix){
		assert(cohort==0 || cohort==1):"Only the two explicitly selected annotation cohorts produce construction FASTAs";
		final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(prefix+suffix,FileFormat.FASTA,null,true,false,false,true),null,1);writer.start();
		ArrayList<Read> batch=new ArrayList<Read>(256);long batchId=0,n=0;
		try{for(Locus l:loci){if(l.cohort!=cohort){continue;}final String id=l.id+" unique="+l.uniqueId+" tid|"+l.taxid+" partial="+l.partial+" pseudo="+l.pseudo;
			batch.add(new Read(l.unique.bases,null,id,n++,false));if(batch.size()==256){writer.add(batch,batchId++);batch=new ArrayList<Read>(256);}
		}if(!batch.isEmpty()){writer.add(batch,batchId);}}finally{require(!writer.poisonAndWait(),"Construction FASTA output failed");}
		require(n==cohortCounts[cohort] && writer.readsWritten()==n,"Export must include every selected copy, including flagged members");
	}
	private static void appendChoice(ByteBuilder b,Choice choice,int preference,int length){
		final int winner=choice.winner();if(winner<0){for(int i=0;i<10;i++){b.append("\tNA");}return;}
		final Hit h=choice.best[winner];b.tab().append(h.model).tab().append(h.bits,9).tab().append(h.coverage,9).tab().append(h.trunc).tab().append(choice.ties()).tab().append(preference<0?"NA":CONCORDANCE[choice.concordance(preference)]);
		if(length>0){final int overlap=choice.overlap[winner];b.tab().append(overlap).tab().append((double)overlap/length,9).tab().append((double)overlap/h.length,9);}else{b.append("\tNA\tNA\tNA");}b.tab().append(h.physical);
	}
	/** Exact token 5S, not 5.8S/15S/25S; quotes and ordinary delimiters are allowed. */
	public static boolean label5s(String product){
		require(product!=null,"A descriptive label must come from the retained source product");
		for(int i=0;i+1<product.length();i++){if(product.charAt(i)!='5' || Character.toLowerCase(product.charAt(i+1))!='s'){continue;}
			if((i==0 || !labelChar(product.charAt(i-1))) && (i+2==product.length() || !labelChar(product.charAt(i+2)))){return true;}}
		return false;
	}
	private static boolean labelChar(char c){return Character.isLetterOrDigit(c) || c=='.';}
	public static int unsupportedCategory(String product,int length,int threshold,boolean anyHit){
		require(length>0 && threshold>0,"Descriptive length categories require measured positive lengths, never an admission cutoff");
		return label5s(product)?0:length<threshold?1:!anyHit?2:3;
	}
	private static int modelIndex(String id){
		for(int i=0;i<MODELS.length;i++){if(MODELS[i].equals(id)){return i;}}
		require(id.equals("RF02540") || id.equals("RF02541") || id.equals("RF02543") || id.equals("RF02546"),"Unexpected CM outside the accepted SSU/LSU roster: "+id);return -1;
	}
	private static int physicalLength(String segments){
		long n=0;for(NcrnaMitoCmEvidence.Span s:NcrnaMitoCmEvidence.parseSegments(segments,Integer.MAX_VALUE)){n+=(long)s.hi-s.lo+1;}
		require(n>0 && n<=Integer.MAX_VALUE,"Physical length must fit one finite evidence record");return (int)n;
	}
	private static boolean bit(LineParser1 p,int field){require(p.termEquals("true",field) || p.termEquals("false",field),"Evidence flags must remain literal booleans");return p.termEquals("true",field);}
	private static void read(String path,int columns,String first,java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{final byte[] header=in.nextLine();require(header!=null,"Missing evidence header: "+path);p.set(header);require(p.terms()==columns && p.termEquals(first,0),"Unexpected accepted evidence schema: "+path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==columns,"Changed evidence row width: "+path);action.accept(p);}}
		finally{require(!in.close(),"Evidence read failed: "+path);}
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Choice {
		void offer(Hit h,int bases){require(h.ssu>=0 && bases>=0,"Only explicit SSU evidence enters the highest-raw-bits comparison");final Hit old=best[h.ssu];
			if(old==null || h.bits>old.bits || h.bits==old.bits && (h.coverage>old.coverage || h.coverage==old.coverage && h.group<old.group)){best[h.ssu]=h;overlap[h.ssu]=bases;}}
		int winner(){int winner=-1;for(int m=0;m<best.length;m++){if(best[m]!=null && (winner<0 || best[m].bits>best[winner].bits)){winner=m;}}return winner;}
		String ties(){final int w=winner();require(w>=0,"Tied-model list needs at least one overlapping SSU hit");final ByteBuilder b=new ByteBuilder();for(int m=0;m<best.length;m++){if(best[m]!=null && best[m].bits==best[w].bits){if(b.length()>0){b.append(',');}b.append(MODELS[m]);}}return b.toString();}
		int concordance(int preference){final int w=winner();require(w>=0 && preference>=0 && preference<3,"Concordance is defined only for annotated universal preferences and an SSU winner");
			if(preference==1){return 4;}int tied=0;for(Hit h:best){if(h!=null && h.bits==best[w].bits){tied++;}}if(tied>1){return 3;}if(w==4){return 2;}return (preference==0)==(w<=1)?0:1;}
		final Hit[] best=new Hit[5];final int[] overlap=new int[5];
	}
	private static final class Hit {int group,ssu,length;long line;String run,model,accession,physical,trunc;double bits,coverage;boolean genomic,crosses,annotated;}
	private static final class Locus {String id,accession,uniqueId,status,product,segments;int taxid,length,preference,cohort;boolean partial,pseudo,crosses,anyHit,lsu,sequenceSeen;Unique unique;final Choice best=new Choice();final double[] identity={Double.NaN,Double.NaN};}
	private static final class Unique {byte[] bases;int copies;final int[] cohortCopies=new int[6];}
	private static final class CmOnly {String accession,physical;int members;final Choice best=new Choice();}
	private String input,locusPath,hitPath,associationPath,pairPath,genomePath,prefix;
	private int expectedLoci=37716,expectedGenomes=18236,expectedHits=89876,expectedGroups=89773,expected16=4831,expected18=6586,expectedUnsupported=8672,shortN;
	private Hit[] hits;private final int[] cohortCounts=new int[6],unsupported=new int[4];
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<CmOnly> cmOnly=new ArrayList<CmOnly>();
	private final ObjectMap<String,Locus> byLocus=new ObjectMap<String,Locus>(String.class,Locus.class);
	private final ObjectMap<String,Unique> byUnique=new ObjectMap<String,Unique>(String.class,Unique.class);
	private final ObjectIntMap<String> taxids=new ObjectIntMap<String>(String.class);
	private static final String[] MODELS={"RF00177","RF01959","RF01960","RF02542","RF02545"};
	private static final String[] PREFERENCES={"16S","TIE","18S"};
	private static final String[] COHORTS={"m16S","m18S","UNRESOLVED_TIE","UNRESOLVED_CONFLICT","UNRESOLVED_UNSUPPORTED","LSU_ONLY"};
	private static final String[] REASONS={"SOURCE_5S_LABEL","SHORTER_THAN_OBSERVED_INTACT_SSU_MINIMUM","NO_GA_HIT_EITHER_STRAND","OPPOSITE_STRAND_ONLY_EVIDENCE"};
	private static final String[] CONCORDANCE={"CONCORDANT_NOMINAL_MODEL_CLASS","DISCORDANT_NOMINAL_MODEL_CLASS","SPECIALIZED_MITO_MODEL","TIED_BEST_MODELS","UNIVERSAL_TIE"};
	private static final String[] SUFFIXES={".audit.tsv.gz",".cross.tsv",".unsupported.tsv",".summary.tsv",".m16S.fa.gz",".m18S.fa.gz"};
}
