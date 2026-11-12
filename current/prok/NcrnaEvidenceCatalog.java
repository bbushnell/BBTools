package prok;

import java.util.ArrayList;
import java.io.File;
import java.nio.charset.StandardCharsets;

import fileIO.FileFormat;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.SequenceInput;
import stream.Read;
import stream.StreamerFactory;
import structures.ListNum;

/** Reads the path-based evaluation census and its original source memberships.
 * Input IDs identify catalog entries, not biological assemblies or TaxIDs.
 * Distinct paths and every source claim survive even when assembly names agree.
 * Loading reads metadata only. Explicit observe calls read the current FASTA;
 * these observations do not prove the bytes used by an earlier scan. Neither
 * operation opens GFF payloads or verifies reported digests. Annotation evidence
 * and biological compartment assignments remain separate from contig lengths.
 * After load returns, catalog maps and memberships are read-only and may be shared
 * by workers. Each observe call owns its reader and mutable Observation; callers
 * must establish native reader globals before starting concurrent observations.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEvidenceCatalog{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads both complete tables, rejecting missing or repeated catalog memberships. */
	public static NcrnaEvidenceCatalog load(final String census, final String memberships){
		if(census==null || memberships==null){
			throw new IllegalArgumentException("Both census and source membership tables are required");
		}
		final NcrnaEvidenceCatalog catalog=new NcrnaEvidenceCatalog();
		catalog.readCensus(census);
		catalog.readMemberships(memberships);
		return catalog;
	}

	private NcrnaEvidenceCatalog(){}

	/*--------------------------------------------------------------*/
	/*----------------           Readers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps the supplied census order; IDs need not be sorted or dense. */
	private void readCensus(final String path){
		assert(genomes.isEmpty()) : "A catalog instance reads one census before its memberships";
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, CENSUS_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==5, "Expected five census fields in "+path);
				final String id=p.parseString(0), fasta=utf8(p, 1);
				final long bytes=p.parseLong(2), bases=p.parseLong(3), records=p.parseLong(4);
				require(validId(id) && fasta.startsWith("/") && bytes>0 && bases>0 && records>0,
					"Census entries need positive IDs/dimensions and absolute source paths: "+id);
				require(!byId.contains(id) && paths.add(fasta),
					"Repeated input ID or resolved FASTA path would double-count evidence: "+id+" "+fasta);
				final Genome genome=new Genome(id, fasta, bytes, bases, records);
				byId.put(id, genome);
				genomes.add(genome);
			}
		}
		require(!genomes.isEmpty(), "An empty census cannot define the evidence denominator: "+path);
	}

	/** Source indices are the original per-source record order, not GFF line numbers. */
	private void readMemberships(final String path){
		assert(!genomes.isEmpty() && membershipCount==0) : "Memberships attach to a complete, previously read census";
		final LineParser1 p=new LineParser1('\t');
		final ObjectIntMap<String> lastIndex=new ObjectIntMap<String>(String.class);
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, MEMBERSHIP_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==10, "Expected ten source membership fields in "+path);
				final String source=utf8(p, 0), id=p.parseString(2);
				final int index=p.parseInt(1);
				final Genome genome=byId.get(id);
				require(!source.isEmpty() && genome!=null,
					"Every source membership must identify a census entry: "+source+" "+id);
				final int previous=lastIndex.contains(source) ? lastIndex.get(source) : 0;
				require(previous<Integer.MAX_VALUE && index==previous+1,
					"Source indices must retain their complete original order: "+source+" "+index);
				final String fasta=utf8(p, 3), sourceGenome=utf8(p, 4), assembly=utf8(p, 5);
				final String taxid=utf8(p, 6), compartment=utf8(p, 7), digest=p.parseString(9);
				p.setBounds(8);
				final String gff=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
				require(fasta.startsWith("/") && !sourceGenome.isEmpty() && !assembly.isEmpty()
					&& !taxid.isEmpty() && !compartment.isEmpty() && !gff.isEmpty(),
					"Missing metadata must remain explicit NA fields, not empty fields: "+source+" "+index);
				require(digest.equals("NA") || SeedOffsetTable.validSha80(digest),
					"A reported digest must remain sha80 or NA; it is not a verified binding: "+source+" "+index);
				genome.members.add(new Membership(source, index, fasta, sourceGenome, assembly, taxid, compartment, gff, digest));
				lastIndex.put(source, index);
				membershipCount++;
			}
		}
		for(final Genome genome:genomes){
			require(!genome.members.isEmpty(), "Census entry lacks source membership: "+genome.id);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Contig observation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads one original FASTA through the native Streamer, retaining only names and lengths.
	 * The caller owns parallelism across inputs. No full-genome sequence copy survives a batch.
	 * Errors in a consumed record are reported after draining the reader so consumer failure
	 * cannot leave its producer blocked on an unread queue. No input hash is manufactured.
	 */
	public Observation observe(final Genome genome){
		require(genome!=null && byId.get(genome.id)==genome, "Observe an entry owned by this catalog");
		require(!Read.FIX_HEADER, "Literal Infernal target identifiers require Read.FIX_HEADER=false before opening FASTA input");
		assert(genome.records>0 && genome.bases>0) : "readCensus supplies positive denominators for every observed input";
		require(new File(genome.fasta).length()==genome.bytes,
			"Source-file size changed since the census: "+genome.fasta);
		final Observation result=new Observation(genome);
		final FileFormat format=FileFormat.testInput(genome.fasta, FileFormat.FASTA, null, true, true);
		require(format.fasta(), "The catalog observation requires FASTA input: "+genome.fasta);
		RuntimeException failure=null;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch.list){
					if(failure!=null){continue;}
					try{result.add(read);}
					catch(RuntimeException invalid){failure=invalid;}
				}
				input.stream.returnList(batch);
			}
			if(failure!=null){throw failure;}
		}
		require(result.size()==genome.records && result.bases==genome.bases,
			"Observed contigs disagree with the census for input "+genome.id+": records="+result.size()
			+"/"+genome.records+", bases="+result.bases+"/"+genome.bases);
		return result;
	}

	/** Starts observations captured by an accepted preparation map, without reopening
	 * the original FASTA. The batch reader must populate every ordinal and validate
	 * totals before publishing this object. The digest is retained provenance, not a
	 * claim that this method recomputed or authenticated the map's bytes.
	 */
	Observation captured(final Genome genome, final String mapping, final String sha80){
		require(genome!=null && byId.get(genome.id)==genome, "Captured observations must belong to this census");
		require(new File(mapping).isAbsolute() && SeedOffsetTable.validSha80(sha80), "Captured observations require their explicit map and sha80");
		return new Observation(genome, mapping, sha80);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of distinct catalog paths, not a count of sequence-unique assemblies. */
	public int size(){return genomes.size();}

	/** Returns one entry in the original census order. */
	public Genome get(final int index){return genomes.get(index);}

	/** Returns the entry for a literal catalog ID, or null if absent. */
	public Genome get(final String id){return byId.get(id);}

	/** Number of retained source rows, including repeated membership of the same input. */
	public long memberships(){return membershipCount;}

	/** Catalog IDs are positive decimal strings with one canonical spelling. */
	private static boolean validId(final String id){
		if(id.isEmpty() || id.charAt(0)<'1' || id.charAt(0)>'9'){return false;}
		for(int i=1; i<id.length(); i++){
			if(id.charAt(i)<'0' || id.charAt(i)>'9'){return false;}
		}
		return true;
	}

	/** Infernal under LC_ALL=C ends identifiers at ASCII whitespace, not arbitrary Unicode spaces. */
	private static String targetToken(final String header){
		require(header!=null && !header.isEmpty(), "FASTA contigs need nonempty identifiers");
		int end=0;
		while(end<header.length()){
			final char c=header.charAt(end);
			if(c==' ' || (c>='\t' && c<='\r')){break;}
			require(c!='\ufffd', "The native reader could not preserve a FASTA identifier's encoding");
			end++;
		}
		require(end>0, "FASTA identifier cannot start with whitespace");
		return end==header.length() ? header : header.substring(0, end);
	}

	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	/** Paths and source claims are retained UTF-8, unlike validated ASCII catalog IDs. */
	private static String utf8(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Metadata            ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal contig observations scoped to one input; equal names in other inputs are unrelated.
	 * Compartment remains a source claim in Membership, not inferred here from a filename or CM hit.
	 */
	public static final class Observation{
		private Observation(final Genome genome_){
			genome=genome_;
			captured=false;
			source=genome.fasta;
			sha80="NA";
		}
		private Observation(final Genome genome_, final String source_, final String sha80_){
			genome=genome_;
			captured=true;
			source=source_;
			sha80=sha80_;
		}

		/** Package-owned capture construction; one canonical ordinal per original record. */
		void addCaptured(final String name, final int length, final long ordinal){
			require(captured && ordinal==names.size()+1L && length>0, "Captured contigs require complete increasing ordinals and positive lengths");
			require(targetToken(name).equals(name) && !lengths.contains(name), "Captured original tokens must be unique within their input");
			names.add(name);
			lengths.put(name, length);
			bases=Math.addExact(bases, length);
		}

		/** No captured observation is complete until it conserves the census dimensions. */
		void finishCaptured(){
			require(captured && size()==genome.records && bases==genome.bases, "Captured map does not cover the exact census contigs/bases for "+genome.id);
		}

		private void add(final Read read){
			assert(names.size()==lengths.size()) : "Each observed literal name has exactly one length before adding another contig";
			require(read!=null && read.bases!=null && read.length()>0 && read.mate==null,
				"Expected one nonempty unpaired FASTA contig in input "+genome.id);
			final String name=targetToken(read.id);
			require(!lengths.contains(name), "Repeated literal contig name in input "+genome.id+": "+name);
			lengths.put(name, read.length());
			names.add(name);
			bases=Math.addExact(bases, read.length());
		}

		public int size(){return names.size();}
		public String name(final int index){return names.get(index);}
		public int length(final int index){return lengths.get(names.get(index));}

		/** Returns -1 for an unobserved literal name; no alias or taxonomy guessing. */
		public int length(final String name){return lengths.contains(name) ? lengths.get(name) : -1;}

		public final Genome genome;
		/** False for a current FASTA read; true for a source-bound preparation map. */
		public final boolean captured;
		/** Actual observation source and its declared map fingerprint (NA for a live read). */
		public final String source, sha80;
		private final ArrayList<String> names=new ArrayList<String>();
		private final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
		private long bases;
	}

	/** One path-based input and its observed whole-genome dimensions. */
	public static final class Genome{
		private Genome(final String id_, final String fasta_, final long bytes_, final long bases_, final long records_){
			id=id_;
			fasta=fasta_;
			bytes=bytes_;
			bases=bases_;
			records=records_;
		}

		public int memberships(){return members.size();}
		public Membership membership(final int index){return members.get(index);}

		public final String id, fasta;
		/** Physical source-file bytes, which may be compressed. */
		public final long bytes;
		/** Observed nucleotide bases and complete FASTA record count. */
		public final long bases, records;
		private final ArrayList<Membership> members=new ArrayList<Membership>();
	}

	/** Original source claims: none is silently promoted to a verified binding. */
	public static final class Membership{
		private Membership(final String source_, final int index_, final String fasta_, final String genome_,
			final String assembly_, final String taxid_, final String compartment_, final String gff_, final String digest_){
			source=source_;
			index=index_;
			suppliedFasta=fasta_;
			sourceGenome=genome_;
			assembly=assembly_;
			taxid=taxid_;
			compartment=compartment_;
			gff=gff_;
			reportedSha80=digest_;
		}

		public final String source;
		public final int index;
		public final String suppliedFasta, sourceGenome, assembly, taxid, compartment;
		/** May be NA or unresolved relative text; reading this record does not open it. */
		public final String gff;
		/** Original reported text; digest kind and actual-file verification remain external evidence. */
		public final String reportedSha80;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ObjectMap<String,Genome> byId=new ObjectMap<String,Genome>(String.class, Genome.class);
	private long membershipCount;

	public static final String CENSUS_HEADER="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records";
	public static final String MEMBERSHIP_HEADER="source\tmember_index\tinput_id\tsupplied_fasta\tsource_genome\tsource_assembly\ttaxid\tdeclared_compartment\treported_gff\treported_sha80";
}
