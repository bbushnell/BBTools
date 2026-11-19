package test.ddl;

import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.lang.reflect.Constructor;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import ddl.DDLFamilyProfiles;
import ddl.DDLKeyTieBreak;
import ddl.FamilySketchBuilder;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import prot.ProteinSequence;
import prot.cluster.ClusterSketchCodec;
import structures.ByteBuilder;

/** Small parent-compatible baseline capture and independent consensus-bonus oracle.
 * Reflection is confined to this fixture so its same source compiles before and
 * after the optional constructor exists; production code uses ordinary calls.
 * @author Brian Bushnell, Sayu
 */
public final class FamilyConsensusBonusFixture{
	public static void main(String[] args) throws Exception{
		boolean assertions=false;
		assert(assertions=true);
		check(assertions && args.length==2,"Use -ea and baseline|candidate plus a fresh output directory");
		final boolean candidate=args[0].equals("candidate");
		check(candidate || args[0].equals("baseline"),"Explicit fixture mode required");
		final Path out=Files.createDirectory(Paths.get(args[1]));
		final byte[] baseline=baselineBytes(false);
		Files.write(out.resolve("default.bin"),baseline);
		if(candidate){
			final byte[] explicit=baselineBytes(true);
			check(Arrays.equals(baseline,explicit),"Explicit W0 must preserve omitted-option profile/query/header bytes");
			Files.write(out.resolve("explicit_zero.bin"),explicit);
			weightedCases();integrationCases(out);
		}
		Files.write(out.resolve("PASS"),("FAMILY_CONSENSUS_BONUS_FIXTURE_PASS mode="+args[0]+" checks="+checks+"\n").getBytes(StandardCharsets.US_ASCII));
		System.out.println("FAMILY_CONSENSUS_BONUS_FIXTURE_PASS mode="+args[0]+" checks="+checks);
	}

	/** Stable binary bytes include old-policy modal choices, query keys and codec headers. */
	private static byte[] baselineBytes(boolean explicitZero) throws Exception{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final DataOutputStream out=new DataOutputStream(bytes);
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(5);
		final byte[][] members={ascii("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"),
			ascii("ACDEFGHIKLMNPQRSAVWYACDEFGHIKLMNPQRSTVWY"),
			ascii("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSAVWY")};
		for(ProteinPositionMode mode:new ProteinPositionMode[]{ProteinPositionMode.PLAIN,ProteinPositionMode.FIXED}){
			for(DDLKeyTieBreak.Policy policy:new DDLKeyTieBreak.Policy[]{DDLKeyTieBreak.Policy.SMALLEST,DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY}){
				final int width=512;
				final DDLKeyTieBreak.PolicyContext context=context(policy,width,mode,0,explicitZero);
				final FamilySketchBuilder builder=new FamilySketchBuilder(width,geometry,0,mode,4,2,null,context);
				for(byte[] member:members){builder.addRaw(member);}
				final byte[] consensus=new ProteinSequence("fixture",members[0]).enc;
				final char[] profile=new char[width];builder.writeProfile(profile,"family.fixed",consensus);
				write(out,profile);
				final FamilySketchBuilder fork=builder.fork();
				for(byte[] member:members){fork.addEncoded(new ProteinSequence("member",member).enc);}
				final char[] rebuilt=new char[width];fork.writeProfile(rebuilt,"family.fixed",consensus);
				check(Arrays.equals(profile,rebuilt),"Fork and raw/encoded member accumulation preserve the same family keys");
				final char[] raw=new char[width],encoded=new char[width];
				builder.sketchRaw(members[0],raw);builder.sketchEncoded(consensus,encoded);
				check(Arrays.equals(raw,encoded),"Raw/encoded query keys agree, including fixed boundary duplication");
				write(out,raw);
				final DDLKeyTieBreak.KeySet all=new DDLKeyTieBreak.KeySet(width);
				DDLFamilyProfiles.enumerateKeysEncoded(consensus,geometry,0,mode,4,2,null,all);
				final char[] maxima=new char[width];all.maxima(maxima);
				check(Arrays.equals(maxima,raw),"Complete consensus enumeration reduces to the ordinary query sketch");
				write(out,maxima);
				final ClusterSketchCodec.Header header=new ClusterSketchCodec.Header(geometry,width,0,mode,4,2,1,null).withKeyPolicy(context);
				final ByteBuilder encodedHeader=new ByteBuilder();ClusterSketchCodec.appendHeader(header,encodedHeader);
				final byte[] headerBytes=encodedHeader.toBytes();out.writeInt(headerBytes.length);out.write(headerBytes);
				builder.reset();check(builder.members()==0 && builder.validWords()==0,"Reset retains no preceding family votes or words");
			}
		}
		for(DDLKeyTieBreak.Policy policy:new DDLKeyTieBreak.Policy[]{DDLKeyTieBreak.Policy.SMALLEST,DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY}){
			final DDLFamilyProfiles.Profile votes=new DDLFamilyProfiles.Profile(4);
			votes.add(new char[]{5,9,0,0});votes.add(new char[]{9,9,0,0});
			final DDLKeyTieBreak.KeySet present=new DDLKeyTieBreak.KeySet(4);
			present.add(0,9);present.add(1,5);present.add(2,7);
			final char[] selected=new char[4];
			votes.writePackedValues(context(policy,4,ProteinPositionMode.PLAIN,0,explicitZero),DDLKeyTieBreak.familySeed("baseline.tie"),present,new DDLKeyTieBreak.Workspace(4),null,selected);
			check(selected[0]==(policy==DDLKeyTieBreak.Policy.SMALLEST?5:9),"Unmodified baseline explicitly exercises consensus versus smaller-key ties");
			check(selected[1]==9 && selected[2]==0,"W0 does not add absent consensus candidates or defeat greater real support");
			write(out,selected);
		}
		out.close();return bytes.toByteArray();
	}

	/** Candidate scores are independently counted, without the production primitive vote map. */
	private static void weightedCases() throws Exception{
		final char[][] members={{9,9,5,0},{9,9,7,0},{9,0,0,0},{7,0,0,0},{7,0,0,0}};
		final int[][] counts=new int[4][16];
		final DDLFamilyProfiles.Profile profile=new DDLFamilyProfiles.Profile(4);
		for(char[] row:members){profile.add(row);for(int b=0;b<4;b++){counts[b][row[b]]++;}}
		final DDLKeyTieBreak.KeySet consensus=new DDLKeyTieBreak.KeySet(4);
		consensus.add(0,7);consensus.add(0,5);consensus.add(1,5);consensus.add(2,5);consensus.add(2,7);
		// Repeated occurrences are presence, not repeated bonus additions.
		consensus.add(0,7);consensus.add(0,5);
		final int[] occupied=new int[4],support=new int[4];
		for(int b=0;b<4;b++){occupied[b]=profile.occupied(b);support[b]=profile.votes(b,false);}
		final long seed=DDLKeyTieBreak.familySeed("bonus.fixed");
		final char[] zero=new char[4];
		profile.writePackedValues(context(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,4,ProteinPositionMode.PLAIN,0,true),seed,consensus,new DDLKeyTieBreak.Workspace(4),null,zero);
		for(long weight:new long[]{1,10,1000,10000}){
			final char[] actual=new char[4];
			final DDLKeyTieBreak.PolicyContext context=context(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,4,ProteinPositionMode.PLAIN,weight,true);
			profile.writePackedValues(context,seed,consensus,new DDLKeyTieBreak.Workspace(4),null,actual);
			boolean changed=false;
			for(int b=0;b<4;b++){
				int best=0,bestTie=-1;long score=-1;
				for(int key=1;key<16;key++){
					final boolean present=consensus.contains(b,key);
					if(counts[b][key]==0 && !present){continue;}
					final long candidate=counts[b][key]+(present?weight:0);
					final int tie=DDLKeyTieBreak.tieSignal(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,key,0,seed,present);
					if(candidate>score || candidate==score && tie>bestTie){best=key;score=candidate;bestTie=tie;}
				}
				check(actual[b]==best,"Independent adjusted-support/tie/canonical-key oracle for W="+weight+" bin="+b);
				check(profile.occupied(b)==occupied[b] && profile.votes(b,false)==support[b],"Bonus does not manufacture real member support/occupancy");
				changed|=actual[b]!=zero[b];
			}
			check(changed && actual[0]==7,"A nonwinning consensus key must change the winning family bin");
			check(actual[1]==(weight==1?9:5),"Consensus-only zero-member-support key enters the candidate set once");
		}
		check(profile.members()==members.length,"Bonus must not add synthetic members");
		boolean rejected=false;
		try{context(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,4,ProteinPositionMode.PLAIN,-1,true);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected,"Negative bonus rejected explicitly");
		rejected=false;
		try{
			profile.writePackedValues(context(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,4,ProteinPositionMode.PLAIN,DDLKeyTieBreak.MAX_VOTES,true),seed,consensus,new DDLKeyTieBreak.Workspace(4),null,new char[4]);
		}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected,"40-bit adjusted-support overflow rejected instead of wrapping");
	}

	/** Positive-W integration, leaving the captured W0 baselineBytes method unchanged. */
	private static void integrationCases(Path out)throws Exception{
		final int width=512;
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(5);
		final byte[] member=new byte[100],consensus=new byte[100];Arrays.fill(member,(byte)'A');Arrays.fill(consensus,(byte)'Y');
		final byte[] encoded=new ProteinSequence("consensus",consensus).enc;
		final DDLKeyTieBreak.PolicyContext weighted=context(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY,width,ProteinPositionMode.FIXED,10,true);
		final FamilySketchBuilder builder=new FamilySketchBuilder(width,geometry,0,ProteinPositionMode.FIXED,4,2,null,weighted);
		builder.addRaw(member);final char[] expected=new char[width];builder.writeProfile(expected,"integration",encoded);
		final FamilySketchBuilder fork=builder.fork();check(fork.policyContext()==weighted && fork.members()==0,"Positive-W fork preserves immutable bonus binding with fresh votes");
		fork.addEncoded(new ProteinSequence("member",member).enc);final char[] values=new char[width];fork.writeProfile(values,"integration",encoded);
		check(Arrays.equals(values,expected),"Positive-W raw/encoded builder fork parity");
		fork.reset();check(fork.members()==0 && fork.validWords()==0,"Positive-W reset clears votes/words");
		fork.addRaw(member);fork.writeProfile(values,"integration",encoded);check(Arrays.equals(values,expected),"Positive-W reset retains no previous-family candidates");
		boolean rejected=false;try{fork.writeProfile(values);}catch(IllegalArgumentException correct){rejected=true;}
		check(rejected,"Positive-W builder requires explicit consensus");
		final ClusterSketchCodec.Header header=new ClusterSketchCodec.Header(geometry,width,0,ProteinPositionMode.FIXED,4,2,1,null).withKeyPolicy(weighted);
		check(header.withFamilies(0).sameKeyPolicy(header),"Roster resizing preserves positive bonus identity");
		final Path codec=out.resolve("positive_bonus.codec");
		ClusterSketchCodec.write(codec.toString(),false,new prot.cluster.Family[0],header.withFamilies(0));
		final ClusterSketchCodec.Header loaded=ClusterSketchCodec.read(codec.toString()).header;
		check(loaded.sameModelBinding(header) && loaded.consensusMemberBonus==10,"Codec write/read retains positive bonus metadata");
		check(!loaded.sameKeyPolicy(new ClusterSketchCodec.Header(geometry,width,0,ProteinPositionMode.FIXED,4,2,0,null)),"W0 and positive-W codecs cannot alias");
		final prot.FamilyMemberPreparation preparation=new prot.FamilyMemberPreparation(geometry,4,2,null,
			prot.FamilyCalibration.MemberRetention.UNIFORM_99,true,prot.MemberLengthRange.LEGACY,
			ProteinPositionMode.FIXED,prot.FamilyAlignmentMetric.LEGACY,width,weighted);
		final prot.FamilyMemberPreparation.Prepared prepared=preparation.prepare(0,"integration",new byte[][]{member},consensus);
		check(Arrays.equals(prepared.positional,expected),"Member preparation uses the same positive-W family selector");
		final char[] memberKeys=new char[width];builder.sketchRaw(member,memberKeys);int zeroSupport=0;
		for(int bin=0;bin<width;bin++){
			final int key=prepared.positional[bin];
			check(prepared.positionalSupport[bin]==(key!=0 && memberKeys[bin]==key?1:0),"Selected-key support counts real members, not W or a displaced modal key");
			check(prepared.positionalOccupancy[bin]==(memberKeys[bin]!=0?1:0),"Positive bonus cannot manufacture member occupancy");
			if(key!=0 && prepared.positionalSupport[bin]==0){zeroSupport++;}
		}
		check(zeroSupport>0,"Integration exercises real zero-member-vote consensus selections");
		rejected=false;
		try{new prot.FamilyMemberPreparation(geometry,4,2,null,prot.FamilyCalibration.MemberRetention.UNIFORM_99,true,
			prot.MemberLengthRange.LEGACY,ProteinPositionMode.PLAIN,prot.FamilyAlignmentMetric.LEGACY,width,
			context(DDLKeyTieBreak.Policy.SMALLEST,width,ProteinPositionMode.PLAIN,10,true));}
		catch(IllegalArgumentException correct){rejected=true;}
		check(rejected,"Unsupported plain positive-W preparation fails before mismatched calibration");
		final Path input=out.resolve("cake_input.faa");Files.write(input,(">member\n"+new String(member,StandardCharsets.US_ASCII)+"\n").getBytes(StandardCharsets.US_ASCII));
		final String[] base={"mode=cake","in="+input,"codecs="+codec,"outdir="+out.resolve("unused_cake"),"fractions=none","matrix=blosum62","memberranking=f"};
		final String[][] options={{},{"consensusbonus=0"},{"consensusbonus=10","consensusbonus=0"},{"consensusbonus=0","consensusbonus=10"}};
		java.util.Map<?,?> omitted=null;
		for(int i=0;i<options.length;i++){
			final String[] args=Arrays.copyOf(base,base.length+options[i].length);System.arraycopy(options[i],0,args,base.length,options[i].length);
			final prot.cluster.DdlClustCake cake=new prot.cluster.DdlClustCake(args,System.err);
			// Test-only inspection avoids adding a production accessor solely for parser regression checks.
			final java.lang.reflect.Field field=prot.cluster.DdlClustCake.class.getDeclaredField("common");field.setAccessible(true);
			final java.util.Map<?,?> actual=(java.util.Map<?,?>)field.get(cake);
			if(i==0){omitted=new java.util.HashMap<Object,Object>(actual);}
			if(i<3){check(actual.equals(omitted) && !actual.containsKey("consensusbonus"),"Omitted, W0 and later-zero Cake options have identical forwarded defaults");}
			else{check("10".equals(actual.get("consensusbonus")),"Last positive Cake option wins");}
		}
	}

	private static DDLKeyTieBreak.PolicyContext context(DDLKeyTieBreak.Policy policy,int width,ProteinPositionMode mode,long bonus,boolean explicit) throws Exception{
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(5);
		if(!explicit){return new DDLKeyTieBreak.PolicyContext(policy,null,width,geometry,0,mode,4,2,null,null);}
		final Constructor<DDLKeyTieBreak.PolicyContext> ctor=DDLKeyTieBreak.PolicyContext.class.getConstructor(
			DDLKeyTieBreak.Policy.class,DDLKeyTieBreak.RawTailCounts.class,int.class,ProteinKmerGeometry.class,
			long.class,ProteinPositionMode.class,int.class,int.class,String.class,String.class,long.class);
		try{return ctor.newInstance(policy,null,width,geometry,0L,mode,4,2,null,null,bonus);}
		catch(java.lang.reflect.InvocationTargetException caught){
			if(caught.getCause() instanceof IllegalArgumentException){throw (IllegalArgumentException)caught.getCause();}
			throw caught;
		}
	}
	private static byte[] ascii(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void write(DataOutputStream out,char[] values) throws Exception{out.writeInt(values.length);for(char v:values){out.writeChar(v);}}
	private static void check(boolean ok,String reason){if(!ok){throw new AssertionError(reason);}checks++;}
	private static int checks;
	private FamilyConsensusBonusFixture(){}
}
