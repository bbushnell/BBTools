package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.LongIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.LongList;

/** Raw heldout coverage by distinct forward17-mer keys; no denominator edits.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSeedCoverage {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaSeedCoverage.class, false);
		String in=null, seeds=null, split=null, out=null, metadata=null, fold="HELD_OUT";int expected=-1;boolean representatives=false, keySupport=false;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}else if(key.equals("seeds")){seeds=value;}else if(key.equals("split")){split=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}
			else if(key.equals("metadata")){metadata=value;}else if(key.equals("fold")){fold=value.toUpperCase(Locale.ROOT);}
			else if(key.equals("representatives")){representatives=Parse.parseBoolean(value);}
			else if(key.equals("keysupport")){keySupport=Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(in!=null && seeds!=null && split!=null && out!=null && expected>0 && (fold.equals("TRAIN") || fold.equals("HELD_OUT") || fold.equals("ALL")),
			"Bind exact selected-fold FASTA, seeds, split ledger, fresh outputs and denominator; fold is TRAIN, HELD_OUT or ALL");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		new NcrnaSeedCoverage().run(in, seeds, split, out, expected, fold, representatives, metadata, keySupport);Shared.closeStream(pp.outstream);
	}

	private void run(final String fasta, final String seeds, final String split, final String out, final int expected, final String fold, final boolean representatives, final String metadata, final boolean keySupport){
		final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
		final ObjectMap<String,Stats> byStratum=new ObjectMap<String,Stats>(String.class, Stats.class);
		final ArrayList<Stats> strata=new ArrayList<Stats>();
		final ByteFile in=RrnaResourceIO.open(split);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group", split);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==7, "Split ledger schema changed");
				require(p.termEquals("TRAIN", 5) || p.termEquals("HELD_OUT", 5), "Unknown split fold");
				if(!fold.equals("ALL") && !p.termEquals(fold, 5) || representatives && !p.termEquals(p.parseString(0), 4)){continue;}
				final String id=p.parseString(0), key=p.parseString(3);require(!members.contains(id), "Duplicate heldout ID: "+id);
				Stats s=byStratum.get(key);if(s==null){s=new Stats(key);byStratum.put(key, s);strata.add(s);}
				members.put(id, new Member(id, p.parseString(2), s));
			}
		}finally{require(!in.close(), "Split ledger read failed");}
		require(members.size()==expected, "Selected-fold population changed: "+members.size()+" expected="+expected);
		final ArrayList<Stats> familyStats=metadata==null?null:loadLabels(metadata, members, expected);
		Collections.sort(strata, (a, b)->a.name.compareTo(b.name));
		final Euk18sSeedMinimumCoverage counter=new Euk18sSeedMinimumCoverage(ProkObject.loadLongKmers(seeds, 17));
		final LongIntMap support=keySupport?new LongIntMap(counter.seeds.size()):null;
		final Stats total=new Stats("ALL");final LongList histogram=new LongList();
		RrnaResourceIO.fresh(out, ".loci.tsv", ".summary.tsv", ".histogram.tsv", ".missed.ids.txt");
		if(metadata!=null){RrnaResourceIO.fresh(out, ".families.tsv");}
		if(keySupport){RrnaResourceIO.fresh(out, ".keys.tsv", ".key_summary.tsv");}
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".loci.tsv"), missed=RrnaResourceIO.writer(out+".missed.ids.txt");
		final ByteBuilder b=new ByteBuilder();
		try{
			rows.println("locus_id\tsource_group\tstratum\tlength\tdistinct_keys");
			RrnaResourceIO.stream(fasta, r->{
				final String id=RrnaResourceIO.first(r.id);final Member m=members.get(id);
				require(m!=null && !m.seen, "FASTA must contain each selected-fold locus exactly once: "+id);m.seen=true;
				final int hits=counter.count(r.bases, support);m.stats.add(hits);total.add(hits);histogram.increment(hits, 1);
				if(m.labels!=null){for(Stats s:m.labels){s.add(hits);}}
				rows.print(b.clear().append(id).tab().append(m.group).tab().append(m.stats.name).tab().append(r.length()).tab().append(hits).nl());
				if(hits==0){missed.println(id);}
			});
		}finally{RrnaResourceIO.finish(rows, missed);}
		require(total.loci==expected, "Missing FASTA members cannot silently shrink the selected-fold denominator: "+total.loci+" expected="+expected);
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv"), hist=RrnaResourceIO.writer(out+".histogram.tsv");
		try{
			summary.println("stratum\tloci\tcovered_ge1\tcovered_ge2\tmissed\tcoverage_ge1\tcoverage_ge2");
			strata.add(0, total);for(Stats s:strata){
				summary.print(b.clear().append(s.name).tab().append(s.loci).tab().append(s.one).tab().append(s.two).tab().append(s.loci-s.one)
					.tab().append(s.one/(double)s.loci, 8).tab().append(s.two/(double)s.loci, 8).nl());
			}
			hist.println("distinct_keys\tloci");for(int i=0; i<histogram.size; i++){
				if(histogram.get(i)>0){hist.print(b.clear().append(i).tab().append(histogram.get(i)).nl());}
			}
		}finally{RrnaResourceIO.finish(summary, hist);}
		if(familyStats!=null){
			final ByteStreamWriter familyOut=RrnaResourceIO.writer(out+".families.tsv");
			try{familyOut.println("view\tfamily\tloci\tcovered_ge1\tcovered_ge2\tmissed\tcoverage_ge1\tcoverage_ge2");for(Stats s:familyStats){
				familyOut.print(b.clear().append(s.name).tab().append(s.loci).tab().append(s.one).tab().append(s.two).tab().append(s.loci-s.one)
					.tab().append(s.one/(double)s.loci, 8).tab().append(s.two/(double)s.loci, 8).nl());
			}}finally{RrnaResourceIO.finish(familyOut);}
		}
		if(keySupport){writeKeySupport(out, counter, support, total, fold, representatives);}
		System.err.println("NCRNA_SEED_COVERAGE_PASS raw="+total.loci+" covered="+total.one+" missed="+(total.loci-total.one));
	}
	private static void writeKeySupport(String out, Euk18sSeedMinimumCoverage counter, LongIntMap support, Stats total, String fold, boolean representatives){
		final long[] keys=counter.seeds.toArray();java.util.Arrays.sort(keys);
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".keys.tsv"), summary=RrnaResourceIO.writer(out+".key_summary.tsv");
		final ByteBuilder b=new ByteBuilder();int zero=0, one=0;long incidences=0;
		try{
			rows.println("forward_17mer\tdistinct_members\tselected_fold_members\tfold\trepresentatives_only");
			for(long key:keys){
				final int count=Math.max(0, support.get(key));
				require(count<=total.loci, "A selected key may contribute at most once per declared input member");
				if(count==0){zero++;}else if(count==1){one++;}incidences+=count;
				rows.print(b.clear().appendKmer(key, 17).tab().append(count).tab().append(total.loci).tab().append(fold).tab().append(representatives).nl());
			}
			require(incidences==total.incidences, "Per-key distinct-member totals must equal summed per-member distinct-key counts");
			summary.println("fold\trepresentatives_only\tmembers\tloaded_keys\tkeys_zero_members\tkeys_one_member\tkeys_below_two_members\tfraction_keys_below_two");
			summary.print(b.clear().append(fold).tab().append(representatives).tab().append(total.loci).tab().append(keys.length).tab().append(zero).tab().append(one)
				.tab().append(zero+one).tab().append((zero+one)/(double)keys.length, 8).nl());
		}finally{RrnaResourceIO.finish(rows, summary);}
	}
	private static ArrayList<Stats> loadLabels(String path, ObjectMap<String,Member> members, int expected){
		final ObjectMap<String,Stats> byLabel=new ObjectMap<String,Stats>(String.class, Stats.class);final ArrayList<Stats> result=new ArrayList<Stats>();
		final String[] views={"SOURCE", "BIOLOGICAL", "BUILD_GROUP"};final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int loaded=0;
		try{RrnaResourceIO.header(in, NcrnaClanSeedPrep.SEED_HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==12, "Family coverage requires complete corrected seed metadata");final Member m=members.get(p.parseString(0));if(m==null){continue;}
			require(m.labels==null, "Each selected member needs one unambiguous label record: "+m.id);m.labels=new Stats[3];loaded++;
			for(int i=0; i<3; i++){
				final String label=p.parseString(9+i);require(!label.isEmpty(), "Selected family label cannot be empty: "+m.id);final String key=views[i]+'\t'+label;
				Stats s=byLabel.get(key);if(s==null){s=new Stats(key);byLabel.put(key, s);result.add(s);}m.labels[i]=s;
			}
		}}finally{require(!in.close(), "Family attribution read failed");}
		require(loaded==expected, "Missing attribution must not shrink per-family denominators");Collections.sort(result, (a, b)->a.name.compareTo(b.name));return result;
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(final String i, final String g, final Stats s){id=i;group=g;stats=s;}
		final String id, group;final Stats stats;boolean seen;Stats[] labels;
	}
	private static final class Stats {
		Stats(final String n){name=n;}
		void add(final int hits){assert(hits>=0):"Distinct loaded-key counts cannot be negative";loci++;incidences+=hits;if(hits>0){one++;}if(hits>1){two++;}}
		final String name;int loci, one, two;long incidences;
	}
}
