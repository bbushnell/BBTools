package test.ddl;

import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import ddl.ProteinKmerWeightedCandidates;
import ddl.ProteinKmerWeightedCandidates.Ranking;
import ddl.ProteinKmerWeightedCandidates.Rankings;
import ddl.ProteinKmerWeightedCandidates.Weights;
import ddl.ProteinKmerGeometry;

/** Independent full-sort rational oracle and aggregation-before-selection fixtures.
 * All synthetic words are software fixtures, never biological evaluation inputs.
 * @author Yaoyao
 */
public final class ProteinKmerWeightedCandidatesTest{
	public static void run(Path root)throws Exception{
		check(!Files.exists(root),"Fresh weighted fixture root required");Files.createDirectory(root);
		weights();projection();oracle(root);aggregation(root);malformed(root);k5(root);dimers(root);
		System.out.println("PROTEIN_KMER_WEIGHTED_CANDIDATES_TEST_PASS checks="+checks);
	}
	private static void weights()throws Exception{
		final Weights w=new Weights(10,1000,4,1);
		check(w.score(3,0)>w.score(0,1000),"Population normalization must override raw-count ranking");
		check(w.score(1,0)==w.score(0,400),"Exact rational tie must remain a tie");
		final Weights same=new Weights(10,1000,8,2);
		check(w.negative==same.negative && w.positive==same.positive,"Equivalent weights must share reduced integer coefficients");
		final Weights actual=new Weights(95000000,1593835520L,4,1);
		check(actual.score(95000000,1593835520L)>0,"Accepted population totals must fit exact long ranking");
		reject(()->new Weights(0,1,4,1),"Zero denominator");
		reject(()->new Weights(1,1,-1,1),"Negative weight");
		reject(()->new Weights(Long.MAX_VALUE,Long.MAX_VALUE-1,4,1),"Overflow must fail, never wrap or approximate");
		reject(()->ProteinKmerWeightedCandidates.parseArms("raw:weighted:2,raw:weighted:02"),"Duplicate equivalent arms");
		for(String a:new String[]{"","raw:shadow:0","reduced:mixed:2","raw:weighted:64000001","raw:weighted:2,"}){
			reject(()->ProteinKmerWeightedCandidates.parseArms(a),"Malformed explicit arms");
		}
		check(ProteinKmerWeightedCandidates.parseArms("raw:shadow:2049,reduced:weighted:16384").length==2,"Sizes are ordinary parameters");
	}
	private static void projection()throws Exception{
		final ProteinKmerWeightedCandidates.Projection p=new ProteinKmerWeightedCandidates.Projection(ALPHABET);
		final ProteinKmerGeometry g=ProteinKmerGeometry.make(6,ALPHABET);
		check(p.universe==2985984,"Twelve-symbol K6 universe");
		final Random random=new Random(61811);
		for(int i=0;i<200;i++){
			final String word=word(random.nextInt(64000000));final int raw=pack(word),dense=reduced(word);
			check(p.project(raw)==dense,"Dense projection against independent symbol oracle");
			check(p.packed(dense)==g.word(raw),"Dense order and production packed geometry must agree");
			check(p.project(p.representative(dense))==dense,"Representative FASTA must round-trip");
		}
		reject(()->new ProteinKmerWeightedCandidates.Projection("all:A/C"),"Incomplete alphabet");
		reject(()->new ProteinKmerWeightedCandidates.Projection("none"),"Implicit alphabet");
	}
	private static void oracle(Path root)throws Exception{
		final Path positive=root.resolve("positive.tsv"),negative=root.resolve("negative.tsv");final Random random=new Random(6137);
		final HashMap<Integer,long[]> raw=new HashMap<Integer,long[]>(),reduced=new HashMap<Integer,long[]>();
		final StringBuilder p=new StringBuilder(HEADER),n=new StringBuilder(HEADER);long np=0,nn=0;int rp=0,rn=0;
		for(int i=0;i<521;i++){
			final String word=word(i);final int key=pack(word);final long pc=i%5==0 ? 0 : 1+random.nextInt(89),nc=i%7==0 ? 0 : 1+random.nextInt(7);
			if(pc>0){row(p,key,word,pc);np+=pc;rp++;}if(nc>0){row(n,key,word,nc);nn+=nc;rn++;}
			if(pc+nc>0){raw.put(key,new long[]{nc,pc});final int r=reduced(word);
				long[] counts=reduced.get(r);if(counts==null){counts=new long[2];reduced.put(r,counts);}counts[0]+=nc;counts[1]+=pc;}
		}
		write(positive,p.toString());write(negative,n.toString());
		for(int limit:new int[]{1,17,83}){
			final Rankings ranks=ProteinKmerWeightedCandidates.rankTables(positive,rp,np,negative,rn,nn,limit,limit,ALPHABET,4,1);
			compare(ranks.rawShadow,raw,false,nn,np,ranks.weights);compare(ranks.rawWeighted,raw,true,nn,np,ranks.weights);
			compare(ranks.reducedShadow,reduced,false,nn,np,ranks.weights);compare(ranks.reducedWeighted,reduced,true,nn,np,ranks.weights);
		}
		final Rankings only=ProteinKmerWeightedCandidates.rankTables(positive,rp,np,negative,rn,nn,0,17,ALPHABET,4,1);
		check(only.rawShadow==null && only.rawWeighted==null,"Reduced-only selection does not allocate raw heaps");
	}
	private static void compare(Ranking actual,HashMap<Integer,long[]> counts,boolean weighted,long nn,long np,Weights weights){
		final ArrayList<Integer> keys=new ArrayList<Integer>();for(int key:counts.keySet()){if(weighted || counts.get(key)[0]>0){keys.add(key);}}
		keys.sort((a,b)->{final int c=score(counts.get(b),weighted,nn,np).compareTo(score(counts.get(a),weighted,nn,np));return c==0 ? Integer.compare(a,b) : c;});
		for(int i=0;i<actual.words.length;i++){
			final int key=keys.get(i);final long[] c=counts.get(key);
			check(actual.words[i]==key,"Exact full-sort rational rank oracle");
			check(actual.negative[i]==c[0] && actual.positive[i]==c[1],"Selected word retains both occurrence counts");
			check(actual.scores[i]==(weighted ? weights.score(c[0],c[1]) : c[0]),"Retained score provenance");
		}
	}
	private static BigInteger score(long[] c,boolean weighted,long nn,long np){
		return weighted ? BigInteger.valueOf(c[0]).multiply(BigInteger.valueOf(np)).multiply(BigInteger.valueOf(4))
			.add(BigInteger.valueOf(c[1]).multiply(BigInteger.valueOf(nn))) : BigInteger.valueOf(c[0]);
	}
	private static void aggregation(Path root)throws Exception{
		final Path p=root.resolve("aggregation_positive.tsv"),n=root.resolve("aggregation_negative.tsv");
		write(p,HEADER+"0\tAAAAAA\t1000\n");
		// D and E each lose to C before reduction, but their combined DE class wins afterward.
		write(n,HEADER+"3\tAAAAAD\t6\n4\tAAAAAC\t10\n5\tAAAAAE\t6\n");
		final Rankings ranks=ProteinKmerWeightedCandidates.rankTables(p,1,1000,n,3,22,1,1,ALPHABET,4,1);
		check(ranks.rawShadow.words[0]==4,"Raw top word is C");
		check(ranks.reducedShadow.words[0]==2 && ranks.reducedShadow.negative[0]==12,"Aggregate synonyms BEFORE ranking");
		check(ranks.reducedWeighted.words[0]==2,"Normalized weighting keeps combined DE above positive-only A");
		write(n,HEADER+"0\tAAAAAA\t1\n1\tAAAAAR\t1\n");write(p,HEADER+"0\tAAAAAA\t1\n1\tAAAAAR\t1\n");
		final Rankings tie=ProteinKmerWeightedCandidates.rankTables(p,2,2,n,2,2,1,1,ALPHABET,4,1);
		check(tie.rawWeighted.words[0]==0 && tie.reducedWeighted.words[0]==0,"Equal scores break by packed key ascending");
	}
	private static void malformed(Path root)throws Exception{
		final Path p=root.resolve("bad_positive.tsv"),n=root.resolve("bad_negative.tsv");write(n,HEADER+"0\tAAAAAA\t1\n");
		for(String row:new String[]{"0\tAAAAAX\t1\n","1\tAAAAAA\t1\n","0\tAAAAAA\t0\n","0\tAAAAAA\t1\n0\tAAAAAA\t1\n","1\tAAAAAR\t1\n0\tAAAAAA\t1\n"}){
			write(p,HEADER+row);reject(()->ProteinKmerWeightedCandidates.rankTables(p,1,1,n,1,1,1,0,null,4,1),"Malformed count table must fail");
		}
		write(p,HEADER+"0\tAAAAAA\t1\n");
		reject(()->ProteinKmerWeightedCandidates.rankTables(p,1,2,n,1,1,1,0,null,4,1),"Window conservation");
		reject(()->ProteinKmerWeightedCandidates.rankTables(p,2,2,n,1,1,1,0,null,4,1),"Row conservation");
		reject(()->ProteinKmerWeightedCandidates.rankTables(p,1,1,n,1,1,2,0,null,4,1),"Insufficient distinct keys");
	}
	/** Independent five-symbol oracle checks leading odd residue, both populations and all ranking modes. */
	private static void k5(Path root)throws Exception{
		final String alphabet="all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T";
		final String[] groups=alphabet.substring(4).split("/");
		final ProteinKmerWeightedCandidates.Projection projection=new ProteinKmerWeightedCandidates.Projection(5,alphabet);
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(5,alphabet);
		check(projection.universe==1048576,"16K5 has exactly 2^20 reduced words");
		final Random random=new Random(16511);
		for(int i=0;i<300;i++){
			final char[] residues=new char[5];int dense=0;
			for(int j=0;j<5;j++){
				residues[j]=AA.charAt(random.nextInt(20));int c=0;while(groups[c].indexOf(residues[j])<0){c++;}dense=dense*16+c;
			}
			final int raw=pack(new String(residues));
			check(projection.project(raw)==dense && projection.packed(dense)==geometry.word(raw),"Odd leading residue agrees with literal symbol oracle");
			check(projection.project(projection.representative(dense))==dense,"K5 FASTA representative round trip");
		}
		final Path p=root.resolve("k5_positive.tsv"),n=root.resolve("k5_negative.tsv");
		// F and W each lose to C before aggregation; their FWY class wins after aggregation.
		write(p,HEADER+"0\tAAAAA\t1000\n4\tAAAAC\t10\n13\tAAAAF\t6\n17\tAAAAW\t6\n");
		write(n,HEADER+"0\tAAAAA\t1\n4\tAAAAC\t10\n13\tAAAAF\t6\n17\tAAAAW\t6\n");
		final HashMap<Integer,long[]> raw=new HashMap<Integer,long[]>(),dense=new HashMap<Integer,long[]>();
		raw.put(0,new long[]{1,1000});raw.put(4,new long[]{10,10});raw.put(13,new long[]{6,6});raw.put(17,new long[]{6,6});
		dense.put(0,new long[]{1,1000});dense.put(1,new long[]{10,10});dense.put(4,new long[]{12,12});
		for(int limit:new int[]{1,2,3}){
			final Rankings ranked=ProteinKmerWeightedCandidates.rankTables(p,4,1022,n,4,23,limit,limit,5,alphabet,4,1);
			compare(ranked.rawShadow,raw,false,23,1022,ranked.weights);compare(ranked.rawWeighted,raw,true,23,1022,ranked.weights);
			compare(ranked.reducedShadow,dense,false,23,1022,ranked.weights);compare(ranked.reducedWeighted,dense,true,23,1022,ranked.weights);
			compareReal(ranked.rawReal,raw);compareReal(ranked.reducedReal,dense);
		}
		final ProteinKmerWeightedCandidates.Arm[] arms=ProteinKmerWeightedCandidates.parseArms("raw:real:1024,reduced:shadow:4096",5);
		check(arms[0].name().equals("raw_real_k5_1024") && arms[0].real && arms[1].name().equals("reduced_shadow_k5_4096"),"K5 names and real-only selection are explicit");
		reject(()->ProteinKmerWeightedCandidates.parseArms("raw:real:3200001",5),"K5 universe overflow");
		reject(()->ProteinKmerWeightedCandidates.parseArms("raw:real:1",4),"Unsupported K");
		write(n,HEADER+"0\tAAAAAA\t1\n");
		reject(()->ProteinKmerWeightedCandidates.rankTables(p,4,1022,n,1,1,1,0,5,null,4,1),"K6 input cannot be silently interpreted as K5");
	}
	/** Three disjoint pair slots: literal substring oracle, complete aggregation and nested requested sizes. */
	private static void dimers(Path root)throws Exception{
		final StringBuilder[] classes=new StringBuilder[120];
		for(int i=0;i<classes.length;i++){classes[i]=new StringBuilder();}
		for(int i=0;i<400;i++){
			final StringBuilder group=classes[i%120];if(group.length()>0){group.append(',');}
			group.append(AA.charAt(i/20)).append(AA.charAt(i%20));
		}
		final String[] groups=new String[120];for(int i=0;i<groups.length;i++){groups[i]=classes[i].toString();}
		final String partition=String.join("/",groups),alphabet="dimer:1,3,5:"+partition;
		final ProteinKmerWeightedCandidates.Projection projection=new ProteinKmerWeightedCandidates.Projection(6,alphabet);
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6,alphabet);
		check(projection.universe==1728000 && projection.alphabet.equals(alphabet),"120^3 keys and exact ordered dimer alphabet");
		final Random random=new Random(1206135);
		for(int i=0;i<2000;i++){
			final String w=word(random.nextInt(64000000));final int dense=dimerDense(w,groups),raw=pack(w);
			check(projection.project(raw)==dense,"Dimer dense keys agree with independent substring/class membership");
			check(projection.packed(dense)==geometry.word(raw),"Selected keys equal production dimer coordinates");
			check(projection.project(projection.representative(dense))==dense,"Representative six-residue FASTA round trip");
		}
		for(int key:new int[]{0,119,120,14399,14400,1727999}){
			final String raw=unpack(projection.representative(key),6);
			check(dimerDense(raw,groups)==key,"Boundary representative agrees with literal oracle");
		}
		final Path positive=root.resolve("dimer_positive.tsv"),negative=root.resolve("dimer_negative.tsv");
		write(positive,HEADER+"0\tAAAAAA\t1000\n");
		write(negative,HEADER+"0\tAAAAAA\t6\n1\tAAAAAR\t10\n192\tAAAAQA\t6\n");
		final Rankings merged=ProteinKmerWeightedCandidates.rankTables(positive,1,1000,negative,3,22,1,1,alphabet,4,1);
		check(merged.rawShadow.words[0]==1 && merged.reducedShadow.words[0]==0 && merged.reducedShadow.negative[0]==12,
			"Dimer synonyms must aggregate BEFORE top-list selection");
		final java.util.TreeMap<Integer,String> words=new java.util.TreeMap<Integer,String>();
		while(words.size()<30000){final String w=word(random.nextInt(64000000));words.put(pack(w),w);}
		final HashMap<Integer,long[]> reduced=new HashMap<Integer,long[]>();
		final StringBuilder ps=new StringBuilder(HEADER),ns=new StringBuilder(HEADER);long np=0,nn=0;int index=0;
		for(java.util.Map.Entry<Integer,String> entry:words.entrySet()){
			final int raw=entry.getKey();final String w=entry.getValue();final long pc=1+index%11,nc=1+index%7;index++;
			row(ps,raw,w,pc);row(ns,raw,w,nc);np+=pc;nn+=nc;
			final int dense=dimerDense(w,groups);long[] counts=reduced.get(dense);
			if(counts==null){counts=new long[2];reduced.put(dense,counts);}counts[0]+=nc;counts[1]+=pc;
		}
		write(positive,ps.toString());write(negative,ns.toString());Ranking previous=null;
		for(int size:new int[]{4096,8192,12288,16384}){
			final Rankings ranks=ProteinKmerWeightedCandidates.rankTables(positive,words.size(),np,negative,words.size(),nn,0,size,alphabet,4,1);
			compare(ranks.reducedShadow,reduced,false,nn,np,ranks.weights);
			compare(ranks.reducedWeighted,reduced,true,nn,np,ranks.weights);
			if(previous!=null){check(Arrays.equals(previous.words,Arrays.copyOf(ranks.reducedShadow.words,previous.words.length)),"All four blacklist sizes are exact ranking prefixes");}
			previous=ranks.reducedShadow;
		}
		reject(()->new ProteinKmerWeightedCandidates.Projection(5,alphabet),"K5 cannot use nonoverlapping K6 representatives");
		reject(()->new ProteinKmerWeightedCandidates.Projection(6,"dimer:1,2,3,4,5:"+partition),"Overlapping representatives need a different algorithm");
		reject(()->new ProteinKmerWeightedCandidates.Projection(6,"dimer:1,3:"+partition),"Incomplete residue coverage");
		reject(()->new ProteinKmerWeightedCandidates.Projection(6,alphabet+",AA"),"Duplicate dimer class membership");
		reject(()->new ProteinKmerWeightedCandidates.Projection(6,"dimer:1,3,5:AA/AR"),"Incomplete dimer partition");
	}
	private static int dimerDense(String word,String[] groups){
		int key=0;
		for(int offset=0;offset<6;offset+=2){
			final String pair=word.substring(offset,offset+2);int found=-1;
			for(int i=0;i<groups.length;i++){if((","+groups[i]+",").contains(","+pair+",")){found=i;break;}}
			check(found>=0,"Independent oracle requires complete pair partition");key=key*groups.length+found;
		}
		return key;
	}
	private static String unpack(int raw,int k){
		final char[] letters=new char[k];for(int i=k-1;i>=0;i--){letters[i]=AA.charAt(raw&31);raw>>>=5;}return new String(letters);
	}

	private static void compareReal(Ranking actual,HashMap<Integer,long[]> counts){
		final ArrayList<Integer> keys=new ArrayList<Integer>();for(int key:counts.keySet()){if(counts.get(key)[1]>0){keys.add(key);}}
		keys.sort((a,b)->{final int c=Long.compare(counts.get(b)[1],counts.get(a)[1]);return c==0 ? Integer.compare(a,b) : c;});
		for(int i=0;i<actual.words.length;i++){
			final int key=keys.get(i);final long[] c=counts.get(key);
			check(actual.words[i]==key && actual.scores[i]==c[1] && actual.negative[i]==c[0] && actual.positive[i]==c[1],"Real-only full-sort/count provenance oracle");
		}
	}

	private static int reduced(String word){int key=0;for(int i=0;i<6;i++){int c=0;while(GROUPS[c].indexOf(word.charAt(i))<0){c++;}key=key*12+c;}return key;}
	private static String word(int value){final char[] w=new char[6];for(int i=5;i>=0;i--){w[i]=AA.charAt(value%20);value/=20;}return new String(w);}
	private static int pack(String word){int key=0;for(int i=0;i<word.length();i++){key=(key<<5)|AA.indexOf(word.charAt(i));}return key;}
	private static void row(StringBuilder s,int key,String word,long count){s.append(key).append('\t').append(word).append('\t').append(count).append('\n');}
	private static void write(Path p,String text)throws Exception{Files.write(p,text.getBytes(StandardCharsets.US_ASCII));}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}checks++;}
	private static void reject(Action action,String why)throws Exception{try{action.run();}catch(IllegalArgumentException | ArithmeticException e){checks++;return;}throw new AssertionError(why);}
	private interface Action{void run()throws Exception;}
	private static int checks;
	private static final String AA="ARNDCEQGHILKMFPSTWYV",ALPHABET="all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST";
	private static final String[] GROUPS=ALPHABET.substring(4).split("/");
	private static final String HEADER="packed_word\tword\toccurrences\n";
}
