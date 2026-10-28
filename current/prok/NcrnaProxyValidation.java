package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins sample provenance to CM evidence without treating unresolved hits as truth.
 * Probability and enriched challenge arms are never pooled into an accuracy rate.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaProxyValidation {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaProxyValidation.class, false);
		String sample=null, audit=null, out=null;int probability=-1, challenge=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("sample")){sample=value;}else if(key.equals("audit")){audit=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expectedprobability")){probability=Parse.parseIntKMG(value);}
			else if(key.equals("expectedchallenge")){challenge=Parse.parseIntKMG(value);}else{throw new IllegalArgumentException(key);}
		}
		require(sample!=null && audit!=null && out!=null && probability>0 && challenge>0, "Bind both evidence tables and each arm's explicit denominator");
		Shared.setThreads(1);new NcrnaProxyValidation().run(sample, audit, out, probability, challenge);Shared.closeStream(pp.outstream);
	}

	private void run(final String sample, final String audit, final String prefix, final int probability, final int challenge){
		load(sample, probability, challenge);
		RrnaResourceIO.fresh(prefix, ".joined.tsv.gz", ".summary.tsv", ".phyla.tsv");
		final ByteFile in=RrnaResourceIO.open(audit);final LineParser1 p=new LineParser1('\t');
		final ByteStreamWriter joined=RrnaResourceIO.writer(prefix+".joined.tsv.gz");final ByteBuilder b=new ByteBuilder();int count=0;
		try{
			RrnaResourceIO.header(in, AUDIT_HEADER, audit);joined.println(AUDIT_HEADER+SAMPLE_SUFFIX);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==34, "CM audit schema changed");
				final String id=p.parseString(0);final Member m=members.get(id);
				require(m!=null && !m.seen, "CM evidence has an unknown or duplicate selected locus: "+id);
				final int end=sourceEnd(row);require(end==m.source.length && samePrefix(row, m.source), "Same locus ID has different source metadata in sample and audit: "+id);
				m.seen=true;count++;
				final String key=m.arm+'\t'+m.eligible+'\t'+m.pilot+'\t'+p.parseString(21)+'\t'+p.parseString(33);
				add(global, totals, key, m.group);add(byPhylum, phyla, m.domain+'\t'+m.phylum+'\t'+key, m.group);
				joined.print(b.clear().append(row).tab().append(m.selection).nl());
			}
		}finally{require(!in.close(), "CM evidence read failed");RrnaResourceIO.finish(joined);}
		require(count==probability+challenge, "Every selected locus must have exactly one CM evidence row, including no-hit loci");
		write(prefix+".summary.tsv", "", totals, probability+challenge);
		write(prefix+".phyla.tsv", "domain\tphylum\t", phyla, probability+challenge);
		System.err.println("NCRNA_PROXY_VALIDATION_PASS probability="+probability+" challenge="+challenge+"; CM evidence categories, not biological FP/FN labels");
	}

	private void load(final String path, final int probability, final int challenge){
		assert(members.size()==0):"A validation report binds exactly one selected sample";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int np=0, nc=0;
		try{
			RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER+SAMPLE_SUFFIX, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==29, "Selected-locus schema changed");
				final String id=p.parseString(0), arm=p.parseString(25);
				require(!members.contains(id), "Duplicate sample locus would bias validation: "+id);
				require(p.termEquals("true", 22) || p.termEquals("false", 22), "Pilot exposure must be explicit");
				require(p.termEquals("true", 23) || p.termEquals("false", 23), "Proxy eligibility must be explicit");
				final boolean pilot=p.termEquals("true", 22), eligible=p.termEquals("true", 23);
				if(arm.equals("PROBABILITY")){require(eligible && !pilot, "Probability arm must retain the declared eligible nonpilot population");np++;}
				else{require(arm.equals("CHALLENGE"), "Unknown sample arm: "+arm);nc++;}
				require(p.termEquals("NOT_ASSIGNED", 26), "This proxy report precedes the final train/holdout split");
				final int end=sourceEnd(row);
				members.put(id, new Member(Arrays.copyOf(row, end), Arrays.copyOfRange(row, end+1, row.length),
					p.parseString(21), p.parseString(4), p.parseString(5), arm, pilot, eligible));
			}
		}finally{require(!in.close(), "Sample selection read failed");}
		require(np==probability && nc==challenge, "Each validation arm must match its declared denominator: "+np+"/"+nc);
	}

	private static int sourceEnd(final byte[] row){
		assert(row!=null):"Source metadata comparisons require the actual retained bytes";
		int tabs=0;for(int i=0;i<row.length;i++){if(row[i]=='\t' && ++tabs==21){return i;}}
		throw new IllegalArgumentException("Missing boundary after21 source metadata fields");
	}
	private static boolean samePrefix(final byte[] row, final byte[] source){
		assert(row.length>source.length):"Audit rows append evidence to the original source metadata";
		for(int i=0;i<source.length;i++){if(row[i]!=source[i]){return false;}}return true;
	}
	private static void add(final ObjectMap<String,Count> map, final ArrayList<Count> list, final String key, final String group){
		assert(group!=null && !group.isEmpty()):"Distinct-source counts require the declared source-group identity";
		Count c=map.get(key);if(c==null){c=new Count(key);map.put(key, c);list.add(c);}c.loci++;c.groups.add(group);
	}
	private static void write(final String path, final String prefix, final ArrayList<Count> rows, final int expected){
		Collections.sort(rows, (a,b)->a.key.compareTo(b.key));
		final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();int total=0;
		try{
			out.println(prefix+"arm\tproxy_eligible\tpilot_exposed_group\tcm_state\tgeometry\tloci\tsource_groups");
			for(Count c:rows){total+=c.loci;out.print(b.clear().append(c.key).tab().append(c.loci).tab().append(c.groups.size()).nl());}
			require(total==expected, "Evidence categories must conserve every selected locus in each summary view");
		}finally{RrnaResourceIO.finish(out);}
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(byte[] s, byte[] t, String g, String d, String p, String a, boolean x, boolean e){source=s;selection=t;group=g;domain=d;phylum=p;arm=a;pilot=x;eligible=e;}
		final byte[] source, selection;final String group, domain, phylum, arm;final boolean pilot, eligible;boolean seen;
	}
	private static final class Count {
		Count(String k){key=k;}final String key;int loci;final ObjectSet<String> groups=new ObjectSet<String>(String.class);
	}
	private final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
	private final ObjectMap<String,Count> global=new ObjectMap<String,Count>(String.class, Count.class), byPhylum=new ObjectMap<String,Count>(String.class, Count.class);
	private final ArrayList<Count> totals=new ArrayList<Count>(), phyla=new ArrayList<Count>();
	public static final String SAMPLE_SUFFIX="\tsource_group\tpilot_exposed_group\tproxy_eligible\tlength_bin\tarm\tsplit\tstratum\tlength_boundary_distance";
	public static final String AUDIT_HEADER=NcrnaCmLocusAudit.HEADER+"\tcm_state\traw_hits\toverlapping_sense_hits\tcm_model_from\tcm_model_to\tcm_seq_from\tcm_seq_to\tcm_trunc\tcm_bits\toverlap_bp\tleft_excess\tright_excess\tgeometry";
}
