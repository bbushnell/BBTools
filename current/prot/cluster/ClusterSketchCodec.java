package prot.cluster;

import java.nio.ByteBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import ddl.DDLKeyTieBreak;
import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prot.Blosum62;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Uncalibrated cluster snapshots: five header lines and five lines per family.
 * Metadata is stable ID/name/member count, followed by consensus, uint16 sketch,
 * float32 monomers, and scalar mask plus five float64 scalars with optional
 * retained-member ranking descriptors. A48 stores raw bits.
 * No acceptance gates, thresholds, member IDs or restart-state claims are invented.
 * @author Brian Bushnell, Fischl
 */
public final class ClusterSketchCodec{

	/*--------------------------------------------------------------*/
	/*----------------             Files            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes caller order unchanged; validates the complete roster before opening output.
	 * Family payloads must remain unchanged throughout this single-producer operation.
	 */
	public static void write(final String path, final boolean overwrite, final Family[] families,
			final int k, final int width, final long sketchSeed, final boolean positional,
			final int partitions, final int overlap){
		write(path,overwrite,families,k,width,sketchSeed,positional,partitions,overlap,null);
	}

	/** Binds the exact optional blacklist bytes used for both family and query sketches.
	 * Null declares no blacklist; the caller obtains the pin from its verified loaded blacklist.
	 */
	public static void write(final String path, final boolean overwrite, final Family[] families,
			final int k, final int width, final long sketchSeed, final boolean positional,
			final int partitions, final int overlap, final String blacklistSha80){
		write(path,overwrite,families,k,width,sketchSeed,positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,
			partitions,overlap,blacklistSha80);
	}

	/** Explicit layout identity prevents relative coordinates from being written as legacy fixed tags. */
	public static void write(final String path, final boolean overwrite, final Family[] families,
			final int k, final int width, final long sketchSeed, final ProteinPositionMode positionMode,
			final int partitions, final int overlap, final String blacklistSha80){
		write(path,overwrite,families,ProteinKmerGeometry.exact(k),width,sketchSeed,positionMode,partitions,overlap,blacklistSha80);
	}

	public static void write(final String path, final boolean overwrite, final Family[] families,
			final Header header){
		writeWithHeader(path,overwrite,families,header.withFamilies(families.length));
	}

	public static void write(final String path, final boolean overwrite, final Family[] families,
			final ProteinKmerGeometry geometry, final int width, final long sketchSeed,
			final ProteinPositionMode positionMode, final int partitions, final int overlap,
			final String blacklistSha80){
		require(path!=null && !path.isEmpty() && families!=null, "Missing cluster output or roster");
		final Header header=new Header(geometry,width,sketchSeed,positionMode,partitions,overlap,families.length,blacklistSha80);
		writeWithHeader(path,overwrite,families,header);
	}

	private static void writeWithHeader(final String path, final boolean overwrite, final Family[] families,
			final Header header){
		final Ids seen=new Ids();
		for(Family family : families){
			validateFamily(family,header);
			seen.add(family.familyId,family.id);
		}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false);
		writer.start();
		try{
			final ByteBuilder buffer=new ByteBuilder();
			appendHeader(header,buffer); writer.print(buffer);
			for(Family family : families){
				buffer.clear(); appendFamily(family,header,buffer); writer.print(buffer);
			}
		}finally{
			if(writer.poisonAndWait()){throw new IllegalStateException("Cluster sketch output failed: "+path);}
		}
	}

	/** Loads only exported data, never fabricating a Family member list from its count. */
	public static Loaded read(final String path){
		require(path!=null && !path.isEmpty(), "Missing cluster sketch input");
		final ByteFile input=ByteFile.makeByteFile(path,false);
		try{
			final byte[][] lines=new byte[LINES][];
			for(int i=0; i<LINES; i++){lines[i]=input.nextLine();}
			final LineParser1 parser=new LineParser1('\t');
			final Header header=parseHeader(lines,parser);
			final ArrayList<Entry> entries=new ArrayList<Entry>(Math.min(header.families,1024));
			final Ids seen=new Ids();
			for(int i=0; i<header.families; i++){
				for(int j=0; j<LINES; j++){lines[j]=input.nextLine();}
				final Entry entry=parseFamily(lines,header,parser);
				seen.add(entry.familyId,entry.name); entries.add(entry);
			}
			require(input.nextLine()==null, "Extra records or trailing data after the declared cluster roster");
			return new Loaded(header,entries.toArray(new Entry[entries.size()]));
		}finally{
			if(input.close()){throw new IllegalStateException("Cluster sketch input failed: "+path);}
		}
	}

	/** Validated read-only view of resident models, without serializing or copying payload arrays.
	 * Call after member publishers join; keep models immutable until fixed-seed import completes.
	 * Optional member descriptors are frozen here, including last-nonempty provenance for empty models.
	 */
	public static Loaded resident(final Header geometry,final Family[] families,final boolean positiveOnly){
		require(geometry!=null && families!=null,"Resident seeds require geometry and a complete roster");
		int count=0;
		for(Family family:families){validateFamily(family,geometry);if(!positiveOnly || family.size()>0){count++;}}
		final Entry[] entries=new Entry[count];
		int next=0;
		for(Family family:families){
			if(positiveOnly && family.size()==0){continue;}
			entries[next++]=new Entry(family.familyId,family.id,family.sketchIdentity(),family.size(),
				family.consensusUnsafe(),family.sketchUnsafe(),family.monomersUnsafe(),family.scalarsUnsafe(),family.scalarMask(),
				family.hasMemberRanking() ? family.memberRankingSnapshot() : null);
		}
		return resident(geometry,entries);
	}

	/** Keeps validated entries borrowed read-only; only the small roster array is owned here. */
	static Loaded resident(final Header geometry,final Entry[] entries){
		require(geometry!=null && entries!=null,"Resident seeds require complete entries and geometry");
		final Ids seen=new Ids();
		for(Entry entry:entries){
			require(entry!=null && entry.familyId>=0 && entry.members>=0,"Resident entry has invalid identity/count");
			validateName(entry.name);
			validateSketchIdentity(entry.sketchIdentity,geometry);
			validatePayload(entry.consensus,entry.sketch,entry.monomers,entry.scalars,entry.scalarMask,geometry);
			validateMemberRanking(entry.members,entry.memberRanking);
			seen.add(entry.familyId,entry.name);
		}
		final Header header=geometry.withFamilies(entries.length);
		return new Loaded(header,entries);
	}

	/** Changes only selected identity; source payload bits/count remain borrowed and unchanged. */
	public static Entry selectedEntry(final Entry source,final int id,final String name){
		require(source!=null && id>=0,"Selected model requires a source and nonnegative global ID");
		validateName(name);
		return new Entry(id,name,source.sketchIdentity,source.members,source.consensus,source.sketch,source.monomers,
			source.scalars,source.scalarMask,source.memberRanking);
	}

	/** Attaches explicitly hydrated original-member statistics without changing identity or model bytes.
	 * The caller binds the source family before renumbering; a count alone cannot prove that binding.
	 */
	public static Entry withMemberRanking(final Entry source,final MemberRankStatistics.Snapshot snapshot){
		require(source!=null && snapshot!=null,"Seed ranking attachment requires a validated source and measured member descriptors");
		return new Entry(source.familyId,source.name,source.sketchIdentity,source.members,source.consensus,source.sketch,source.monomers,
			source.scalars,source.scalarMask,snapshot);
	}

	/** Replaces only a measured sketch, retaining all original identity/count/model fields.
	 * Payload arrays are borrowed read-only, as with selectedEntry; the caller must bind
	 * original members to this entry and use the declared geometry to build the new sketch.
	 */
	public static Entry withSketch(final Entry source, final char[] sketch, final Header header){
		require(source!=null && header!=null, "Sketch replacement requires a source entry and target geometry");
		require(header.keyPolicy==DDLKeyTieBreak.Policy.SMALLEST && header.consensusMemberBonus==0,
			"Nondefault key-policy sketches need an explicit effective sketch identity");
		return withSketch(source,sketch,source.name,header);
	}

	public static Entry withSketch(final Entry source, final char[] sketch, final String sketchIdentity,
			final Header header){
		require(source!=null && header!=null, "Sketch replacement requires a source entry and target geometry");
		validateSketchIdentity(sketchIdentity,header);
		validatePayload(source.consensus, sketch, source.monomers, source.scalars, source.scalarMask, header);
		return new Entry(source.familyId, source.name, sketchIdentity, source.members, source.consensus, sketch,
			source.monomers, source.scalars, source.scalarMask, source.memberRanking);
	}

	/** Applies one structural edit without copying payload arrays or serializing the full base.
	 * Loaded rosters originate from read(), which already validates their unique IDs/names.
	 * Neither input roster nor any Entry/payload is modified. Survivors keep base order.
	 */
	public static Loaded applyOverlay(final Loaded base,final Loaded overlay,final int[] removed){
		require(base!=null && overlay!=null && removed!=null && (removed.length==1 || removed.length==2),
			"Seed overlay requires a base, replacements and one or two explicit removed IDs");
		final Header a=base.header,b=overlay.header;
		require(a.families==base.families.length && b.families==overlay.families.length,
			"Seed overlay roster census differs from its header");
		require(a.sameModelBinding(b), "Seed overlay model binding differs from the original seeds");
		for(int i=0;i<removed.length;i++){
			require(removed[i]>=0 && removed[i]<Integer.MAX_VALUE,"Removed seed IDs must be nonnegative and bounded");
			for(int j=0;j<i;j++){require(removed[i]!=removed[j],"Duplicate removed seed ID");}
		}
		final boolean fission=removed.length==1;
		require(overlay.families.length==(fission ? 2 : 1),"Seed overlay must describe exactly one fission or fusion");
		final int retainedId=fission ? removed[0] : Math.min(removed[0],removed[1]);
		final Ids seen=new Ids();
		Entry retained=null,fresh=null;
		long replacementMembers=0,removedMembers=0;
		for(Entry entry:overlay.families){
			require(entry!=null && entry.familyId>=0 && entry.familyId<Integer.MAX_VALUE && entry.members>0,
				"Replacement seeds need bounded IDs and positive member counts");
			seen.add(entry.familyId,entry.name);
			replacementMembers=Math.addExact(replacementMembers,(long)entry.members);
			if(entry.familyId==retainedId){retained=entry;}else{fresh=entry;}
		}
		require(retained!=null,"Seed overlay must retain the parent ID or the smaller fusion parent ID");
		final int[] found=new int[removed.length];
		for(Entry entry:base.families){
			require(entry!=null,"Missing original seed entry");
			boolean remove=false;
			for(int j=0;j<removed.length;j++){
				if(entry.familyId==removed[j]){
					require(++found[j]==1 && entry.members>0,"Removed seed ID must exist exactly once with positive members");
					removedMembers=Math.addExact(removedMembers,(long)entry.members); remove=true;
				}
			}
			if(fresh!=null){require(entry.familyId!=fresh.familyId,"Fission right-child ID must be fresh in the original roster");}
			if(!remove){
				for(Entry replacement:overlay.families){
					require(entry.familyId!=replacement.familyId && !entry.name.equals(replacement.name),
						"Seed overlay collides with a surviving seed ID or name");
				}
			}
		}
		for(int count:found){require(count==1,"Removed seed ID is absent from the original roster");}
		require(removedMembers==replacementMembers,"Seed overlay must conserve removed original-member counts");
		final Entry[] entries=new Entry[Math.addExact(base.families.length-removed.length,overlay.families.length)];
		int next=0;
		for(Entry entry:base.families){
			if(entry.familyId==retainedId){entries[next++]=retained;}
			else if(removed.length==1 || entry.familyId!=removed[0] && entry.familyId!=removed[1]){entries[next++]=entry;}
		}
		if(fresh!=null){entries[next++]=fresh;}
		require(next==entries.length,"Seed overlay output census mismatch");
		return new Loaded(a.withFamilies(entries.length),entries);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Encoding           ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends one fixed five-line header, including explicit absence of calibration. */
	public static void appendHeader(final Header h, final ByteBuilder out){
		require(h!=null && out!=null, "Missing cluster header/output");
		out.append("#cluster_sketch").tab().append(SCHEMA).tab().append("uncalibrated").tab()
			.append(h.k).tab().append(h.width).tab().append(h.sketchSeed).tab()
			.append(h.positionMode.sketch).tab()
			.append(h.partitions).tab().append(h.overlap).tab().append(h.families).tab()
			.append(h.blacklistSha80==null ? "." : h.blacklistSha80);
		if(!h.geometry.isExact()){out.tab().append(h.geometry.name());}
		if(h.keyPolicy!=DDLKeyTieBreak.Policy.SMALLEST || h.consensusMemberBonus>0){
			out.tab().append("keypolicy=").append(h.keyPolicy.label);
			out.tab().append("keyencoding=").append(h.keyEncodingVersion);
			if(h.keyRawInputBinding!=null){out.tab().append("keyraw="); appendEscaped(h.keyRawInputBinding,out);}
		}
		if(h.consensusMemberBonus>0){out.tab().append("consensusbonus=").append(h.consensusMemberBonus);}
		out.nl();
		for(byte[] label : LABELS){out.append(label).nl();}
	}

	/** Appends five family lines; numeric arrays retain their exact IEEE and unsigned bits. */
	public static void appendFamily(final Family family, final Header header, final ByteBuilder out){
		require(out!=null, "Missing cluster record output");
		validateFamily(family,header);
		out.append("#family").tab().append(family.familyId).tab().append(family.id).tab().append(family.size());
		if(header.keyPolicy!=DDLKeyTieBreak.Policy.SMALLEST || header.consensusMemberBonus>0){
			out.tab().append("sketchid="); appendEscaped(family.sketchIdentity(),out);
		}
		out.nl();
		appendPayload(family.consensusUnsafe(),family.sketchUnsafe(),family.monomersUnsafe(),
			family.scalarsUnsafe(),family.scalarMask(),
			family.hasMemberRanking() ? family.memberRankingSnapshot() : null,out);
	}

	/** Writes a selected decoded model with a new stable identity but unchanged payload/count.
	 * No Family or fictional membership list is constructed from the source count.
	 */
	public static void appendEntry(final Entry entry,final int familyId,final String name,
			final Header header,final ByteBuilder out){
		require(entry!=null && header!=null && out!=null && familyId>=0 && entry.members>=0,
			"Selected entry needs valid metadata, header and output");
		validateName(name);
		validateSketchIdentity(entry.sketchIdentity,header);
		validatePayload(entry.consensus,entry.sketch,entry.monomers,entry.scalars,entry.scalarMask,header);
		out.append("#family").tab().append(familyId).tab().append(name).tab().append(entry.members);
		if(header.keyPolicy!=DDLKeyTieBreak.Policy.SMALLEST || header.consensusMemberBonus>0){
			out.tab().append("sketchid="); appendEscaped(entry.sketchIdentity,out);
		}
		out.nl();
		appendPayload(entry.consensus,entry.sketch,entry.monomers,entry.scalars,entry.scalarMask,entry.memberRanking,out);
	}

	/** Shared encoder preserves exact float and scalar bits for rebuilt and selected models. */
	private static void appendPayload(final byte[] consensus,final char[] sketch,final float[] monomers,
			final double[] scalars,final int scalarMask,
			final MemberRankStatistics.Snapshot memberRanking,final ByteBuilder out){
		assert(out!=null && consensus!=null && sketch!=null && monomers!=null && scalars!=null) :
			"Payload arrays must pass codec validation before encoding";
		for(byte residue : consensus){
			out.append(residue==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[residue]);
		}
		out.nl();
		for(int i=0; i<sketch.length; i++){if(i>0){out.tab();} out.appendA48((long)sketch[i]);}
		out.nl();
		for(int i=0; i<monomers.length; i++){if(i>0){out.tab();} out.appendFloatA48(monomers[i]);}
		out.nl().appendA48(scalarMask);
		for(double scalar : scalars){out.tab().appendA48(Double.doubleToRawLongBits(scalar));}
		if(memberRanking!=null){
			out.tab().append(memberRanking.members).tab()
				.appendA48(Double.doubleToRawLongBits(memberRanking.medianLength)).tab()
				.appendA48(Double.doubleToRawLongBits(memberRanking.meanKd)).tab()
				.append(memberRanking.layer);
		}
		out.nl();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Decoding           ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires exact schema/labels and validates geometry before consuming payload rows. */
	public static Header parseHeader(final byte[][] lines, final LineParser1 parser){
		requireLines(lines); require(parser!=null, "Missing header parser");
		parser.set(lines[0]);
		require(parser.terms()>=11 && parser.termEquals("#cluster_sketch",0) && parser.termEquals(SCHEMA,1) &&
			parser.termEquals("uncalibrated",2), "Unsupported cluster sketch schema or calibration declaration");
		final int k=integer(parser,3);
		int next=11;
		ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(k);
		if(parser.terms()>next && !hasEquals(parser,next)){geometry=ProteinKmerGeometry.fromName(parser.parseString(next++));}
		final DDLKeyTieBreak.Policy[] keyPolicy={DDLKeyTieBreak.Policy.SMALLEST};
		final String[] keyEncoding={DDLKeyTieBreak.encodingVersion(keyPolicy[0])}, keyRaw={null};
		final boolean[] seenPolicy={false}, seenEncoding={false}, seenRaw={false};
		long consensusMemberBonus=0;boolean seenBonus=false;
		for(; next<parser.terms(); next++){
			parser.setBounds(next);
			final byte[] line=parser.line();
			int equals=parser.a();
			while(equals<parser.b() && line[equals]!='='){equals++;}
			require(equals>parser.a() && equals<parser.b(), "Malformed cluster key-policy metadata");
			final String value=new String(line,equals+1,parser.b()-equals-1,StandardCharsets.UTF_8);
			if(termEquals(line,parser.a(),equals,"keypolicy")){
				require(!seenPolicy[0], "Duplicate key-policy metadata");
				keyPolicy[0]=DDLKeyTieBreak.Policy.parse(value); seenPolicy[0]=true;
			}else if(termEquals(line,parser.a(),equals,"keyencoding")){
				require(!seenEncoding[0], "Duplicate key-encoding metadata");
				keyEncoding[0]=value; seenEncoding[0]=true;
			}else if(termEquals(line,parser.a(),equals,"keyraw")){
				require(!seenRaw[0], "Duplicate key-raw metadata");
				keyRaw[0]=unescape(value); seenRaw[0]=true;
			}else if(termEquals(line,parser.a(),equals,"consensusbonus")){
				require(!seenBonus,"Duplicate consensus bonus metadata");
				consensusMemberBonus=Long.parseLong(value);seenBonus=true;
			}else{throw new IllegalArgumentException("Unknown cluster key-policy metadata");}
		}
		if(seenPolicy[0] || seenEncoding[0] || seenRaw[0] || seenBonus){
			require(seenPolicy[0] && seenEncoding[0], "Explicit key metadata requires keypolicy and keyencoding");
		}
		require(geometry.k()==k, "Cluster header geometry must match the declared base k");
		final ProteinPositionMode positionMode=ProteinPositionMode.fromSketch(parser.parseString(6));
		final Header header=new Header(geometry,integer(parser,4),signedLong(parser,5),positionMode,
				integer(parser,7),integer(parser,8),integer(parser,9),parser.termEquals(".",10) ? null : parser.parseString(10),
				keyPolicy[0],keyEncoding[0],keyRaw[0],consensusMemberBonus);
		for(int i=0; i<LABELS.length; i++){require(Arrays.equals(lines[i+1],LABELS[i]), "Unknown cluster payload encoding/label");}
		return header;
	}

	/** Parses one complete record into owned data arrays, with no normalization or gate inference. */
	public static Entry parseFamily(final byte[][] lines, final Header header, final LineParser1 parser){
		requireLines(lines); require(header!=null && parser!=null, "Missing cluster header/parser");
		parser.set(lines[0]);
		require((parser.terms()==4 || parser.terms()==5) && parser.termEquals("#family",0),
			"Expected stable ID, name and member count");
		final int familyId=integer(parser,1), members=integer(parser,3);
		final String name=parser.parseString(2); validateName(name);
		final String sketchIdentity;
		if(parser.terms()==5){
			final String metadata=parser.parseString(4);
			require(metadata.startsWith("sketchid="), "Unknown cluster family metadata");
			sketchIdentity=unescape(metadata.substring(9));
		}else{
			require(header.keyPolicy==DDLKeyTieBreak.Policy.SMALLEST && header.consensusMemberBonus==0,
				"Nondefault key-policy records require explicit sketch identity");
			sketchIdentity=name;
		}
		validateSketchIdentity(sketchIdentity,header);
		final byte[] raw=lines[1];
		require(raw.length>0, "A cluster consensus cannot be empty");
		final byte[] consensus=new byte[raw.length];
		for(int i=0; i<raw.length; i++){
			final int residue=raw[i], code=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
			if(code>=0 && code<20 && AminoAcid.numberToAcid[code]==residue){consensus[i]=(byte)code;}
			else if(residue=='X'){consensus[i]=Blosum62.X_CODE;}
			else{throw new IllegalArgumentException("Consensus must contain canonical uppercase standard residues or X");}
		}
		parser.set(lines[2]); require(parser.terms()==header.width, "Sketch width differs from the header");
		final char[] sketch=new char[header.width];
		for(int i=0; i<sketch.length; i++){sketch[i]=(char)a48(parser,i,16);}
		parser.set(lines[3]); require(parser.terms()==MONOMERS, "Expected twenty monomer values");
		final float[] monomers=new float[MONOMERS];
		for(int i=0; i<monomers.length; i++){
			monomers[i]=Float.intBitsToFloat((int)a48(parser,i,32));
			require(Float.isFinite(monomers[i]) && monomers[i]>=0, "Monomers must be finite and nonnegative");
		}
		parser.set(lines[4]);
		final int scalarTerms=parser.terms();
		require(scalarTerms==SCALARS+1 || scalarTerms==SCALARS+1+RANKING_FIELDS,
			"Expected scalar mask, five scalar values and optional member-ranking descriptors");
		final int mask=(int)a48(parser,0,5);
		final double[] scalars=new double[SCALARS];
		for(int i=0; i<scalars.length; i++){
			scalars[i]=Double.longBitsToDouble(a48(parser,i+1,64));
			require(Double.isFinite(scalars[i]), "Scalars must be finite even when their mask bit is clear");
		}
		MemberRankStatistics.Snapshot memberRanking=null;
		if(scalarTerms==SCALARS+1+RANKING_FIELDS){
			memberRanking=new MemberRankStatistics.Snapshot(integer(parser,SCALARS+1),
				Double.longBitsToDouble(a48(parser,SCALARS+2,64)),
				Double.longBitsToDouble(a48(parser,SCALARS+3,64)),integer(parser,SCALARS+4));
		}
		return new Entry(familyId,name,sketchIdentity,members,consensus,sketch,monomers,scalars,mask,memberRanking);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks wire invariants without allocating per-field or per-array conversions. */
	private static void validateFamily(final Family family, final Header header){
		require(family!=null && header!=null, "Missing family/header");
		require(family.familyId>=0 && family.size()>=0, "Negative stable family ID or member count");
		validateName(family.id);
		validateSketchIdentity(family.sketchIdentity(),header);
		validatePayload(family.consensusUnsafe(),family.sketchUnsafe(),family.monomersUnsafe(),
			family.scalarsUnsafe(),family.scalarMask(),header);
	}

	/** Enforces identical wire invariants regardless of the source of the payload. */
	private static void validatePayload(final byte[] consensus,final char[] sketch,final float[] monomers,
			final double[] scalars,final int scalarMask,final Header header){
		require(consensus!=null && consensus.length>0 && sketch!=null && sketch.length==header.width &&
			monomers!=null && monomers.length==MONOMERS && scalars!=null && scalars.length==SCALARS,
			"Family payload dimensions disagree with the cluster sketch schema");
		for(byte residue : consensus){require((residue>=0 && residue<20) || residue==Blosum62.X_CODE, "Invalid encoded consensus residue");}
		for(float value : monomers){require(Float.isFinite(value) && value>=0, "Monomers must be finite and nonnegative");}
		for(double value : scalars){require(Double.isFinite(value), "Scalars must be finite");}
		require(scalarMask>=0 && scalarMask<(1<<SCALARS), "Scalar validity mask exceeds five bits");
	}

	/** Names are retained as-is and cannot inject whitespace into record metadata. */
	private static void validateName(final String name){
		require(name!=null && !name.isEmpty(), "Missing cluster name");
		for(int i=0; i<name.length(); i++){
			final char value=name.charAt(i);
			require(value>32 && value<127, "Cluster names must be printable ASCII tokens");
		}
	}

	private static void validateSketchIdentity(final String identity, final Header header){
		validateName(identity);
	}

	/** Detects optional key/value metadata after the legacy fixed header fields. */
	private static boolean hasEquals(final LineParser1 parser, final int term){
		parser.setBounds(term);
		for(int i=parser.a(); i<parser.b(); i++){if(parser.line()[i]=='='){return true;}}
		return false;
	}

	/** Compares a raw metadata key without allocating a temporary String. */
	private static boolean termEquals(final byte[] line, final int a, final int b, final String expected){
		if(b-a!=expected.length()){return false;}
		for(int i=0; i<expected.length(); i++){if(line[a+i]!=(byte)expected.charAt(i)){return false;}}
		return true;
	}

	/** Percent-escapes arbitrary UTF-8 metadata into one printable tab-safe token. */
	private static void appendEscaped(final String value, final ByteBuilder out){
		require(value!=null && !value.isEmpty(), "Escaped metadata cannot be empty");
		for(byte b : value.getBytes(StandardCharsets.UTF_8)){
			final int c=b&255;
			if(c>32 && c<127 && c!='%'){
				out.append((char)c);
			}else{
				out.append('%').append(HEX[c>>>4]).append(HEX[c&15]);
			}
		}
	}

	/** Reverses percent escapes with strict UTF-8; malformed bytes must not normalize. */
	private static String unescape(final String value){
		require(value!=null && !value.isEmpty(), "Escaped metadata cannot be empty");
		final byte[] out=new byte[value.length()];
		int length=0;
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c!='%'){
				require(c>32 && c<127, "Unescaped metadata must be a printable ASCII token");
				out[length++]=(byte)c;
			}else{
				require(i+2<value.length(), "Truncated percent escape in metadata");
				final int hi=hex(value.charAt(++i)), lo=hex(value.charAt(++i));
				out[length++]=(byte)((hi<<4)|lo);
			}
		}
		try{
			return StandardCharsets.UTF_8.newDecoder().onMalformedInput(CodingErrorAction.REPORT).
				onUnmappableCharacter(CodingErrorAction.REPORT).decode(ByteBuffer.wrap(out,0,length)).toString();
		}catch(CharacterCodingException e){throw new IllegalArgumentException("Malformed UTF-8 in metadata",e);}
	}

	private static int hex(final char c){
		if(c>='0' && c<='9'){return c-'0';}
		if(c>='A' && c<='F'){return c-'A'+10;}
		if(c>='a' && c<='f'){return c-'a'+10;}
		throw new IllegalArgumentException("Invalid percent escape in metadata");
	}

	/** ByteFile.nextLine returns null at EOF; a partial record cannot be mistaken for an empty roster. */
	private static void requireLines(final byte[][] lines){
		require(lines!=null && lines.length==LINES, "Expected a complete five-line block");
		for(byte[] line : lines){require(line!=null, "Truncated cluster header or family record");}
	}

	/** Checks native A48 width before decoding. ByteBuilder.appendA48(long) sign-extends
	 * negative 64-bit patterns into its first 6-bit symbol (56..63, not unsigned 8..15).
	 * Requiring that extension prevents accepting arbitrary overflowing 66-bit tokens.
	 */
	private static long a48(final LineParser1 parser, final int term, final int bits){
		final int length=parser.setBounds(term), start=parser.a(), symbols=(bits+5)/6;
		final byte[] line=parser.line();
		require(length>0 && length<=symbols && (length==1 || line[start]!='0'), "Noncanonical A48 field length/leading zero");
		for(int i=start; i<parser.b(); i++){require(line[i]>=48 && line[i]<=111, "Invalid A48 alphabet");}
		final int first=line[start]-48, topBits=bits-6*(symbols-1);
		if(bits==64 && length==11){
			require(first<8 || first>=56, "A48 float64 pattern has an invalid sign extension");
		}else{
			require(length<symbols || first<(1<<topBits), "A48 field overflows its declared bit width");
		}
		return parser.parseLongA48(term);
	}

	/** Decimal metadata uses overflow-checked accumulation rather than wrapping int parsing. */
	private static int integer(final LineParser1 parser, final int term){
		final long value=signedLong(parser,term);
		require(value>=0 && value<=Integer.MAX_VALUE, "Cluster integer is outside the nonnegative int range");
		return (int)value;
	}

	/** Negative accumulation supports every signed sketch seed, including Long.MIN_VALUE. */
	private static long signedLong(final LineParser1 parser, final int term){
		final int length=parser.setBounds(term); require(length>0, "Empty decimal metadata");
		final byte[] line=parser.line(); int pos=parser.a();
		final boolean negative=line[pos]=='-'; if(negative){pos++;}
		require(pos<parser.b(), "Missing decimal metadata digits");
		long value=0; final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
		for(; pos<parser.b(); pos++){
			final int digit=line[pos]-'0';
			require(digit>=0 && digit<=9 && value>=limit/10, "Invalid or overflowing decimal metadata");
			value*=10; require(value>=limit+digit, "Decimal metadata overflow"); value-=digit;
		}
		return negative ? value : -value;
	}

	private static void require(final boolean ok, final String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	/** Empty templates retain a previous nonempty population; nonempty entries must describe their current members. */
	private static void validateMemberRanking(final int members,final MemberRankStatistics.Snapshot snapshot){
		require(members>=0 && (snapshot==null || members==0 || snapshot.members==members),
			"Resident member descriptors must cover the exact nonempty roster or retain an empty template's last-nonempty population");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Data Types          ----------------*/
	/*--------------------------------------------------------------*/

	/** Version-two geometry, blacklist identity, and optional nondefault key policy. */
	public static final class Header{
		public Header(final int k_, final int width_, final long sketchSeed_, final boolean positional_,
				final int partitions_, final int overlap_, final int families_){
			this(k_,width_,sketchSeed_,positional_,partitions_,overlap_,families_,null);
		}

		/** Null is the only API representation of no blacklist; wire encoding uses a dot. */
		public Header(final int k_, final int width_, final long sketchSeed_, final boolean positional_,
				final int partitions_, final int overlap_, final int families_, final String blacklistSha80_){
			this(k_,width_,sketchSeed_,positional_ ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,
				partitions_,overlap_,families_,blacklistSha80_);
		}

		/** Positional is a compatibility view; positionMode carries the complete stored identity. */
		public Header(final int k_, final int width_, final long sketchSeed_, final ProteinPositionMode positionMode_,
				final int partitions_, final int overlap_, final int families_, final String blacklistSha80_){
			this(ProteinKmerGeometry.exact(k_),width_,sketchSeed_,positionMode_,partitions_,overlap_,families_,blacklistSha80_);
		}

		public Header(final ProteinKmerGeometry geometry_, final int width_, final long sketchSeed_, final ProteinPositionMode positionMode_,
				final int partitions_, final int overlap_, final int families_, final String blacklistSha80_){
			this(geometry_,width_,sketchSeed_,positionMode_,partitions_,overlap_,families_,blacklistSha80_,
				DDLKeyTieBreak.Policy.SMALLEST,DDLKeyTieBreak.encodingVersion(DDLKeyTieBreak.Policy.SMALLEST),null);
		}

		public Header(final ProteinKmerGeometry geometry_, final int width_, final long sketchSeed_, final ProteinPositionMode positionMode_,
				final int partitions_, final int overlap_, final int families_, final String blacklistSha80_,
				final DDLKeyTieBreak.Policy keyPolicy_, final String keyEncodingVersion_, final String keyRawInputBinding_){
			this(geometry_,width_,sketchSeed_,positionMode_,partitions_,overlap_,families_,blacklistSha80_,
				keyPolicy_,keyEncodingVersion_,keyRawInputBinding_,0);
		}

		/** Explicit family-only selection bonus; omitted/default metadata remains byte-identical. */
		public Header(final ProteinKmerGeometry geometry_, final int width_, final long sketchSeed_, final ProteinPositionMode positionMode_,
				final int partitions_, final int overlap_, final int families_, final String blacklistSha80_,
				final DDLKeyTieBreak.Policy keyPolicy_, final String keyEncodingVersion_, final String keyRawInputBinding_,
				final long consensusMemberBonus_){
			require(consensusMemberBonus_>=0 && consensusMemberBonus_<=DDLKeyTieBreak.MAX_VOTES,"Invalid consensus member bonus");
			require(positionMode_!=null, "Cluster headers require an explicit position layout");
			require(geometry_!=null && width_>=1 && width_<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width_)==1,
				"Invalid cluster geometry or sketch width");
			require(partitions_>=0 && overlap_>=0 && families_>=0, "Negative cluster geometry or roster count");
			require(keyPolicy_!=null && DDLKeyTieBreak.encodingVersion(keyPolicy_).equals(keyEncodingVersion_),
				"Cluster header key policy and encoding version differ");
			require((keyPolicy_==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY)==(keyRawInputBinding_!=null),
				"Only rawtail key policy carries a raw input binding");
			if(blacklistSha80_!=null){
				require(geometry_.span()>=5 && (blacklistSha80_.length()==20 || blacklistSha80_.equals("unmeasured")),
					"A filtered codec requires a five-residue span and a sha80 or explicit unmeasured provenance");
				// Unmeasured is never content equivalence. Structural experiment controllers
				// must bind the exact same immutable filter for every seed/continuation.
				for(int i=0; !blacklistSha80_.equals("unmeasured") && i<blacklistSha80_.length(); i++){
					final char c=blacklistSha80_.charAt(i);
					require((c>='0' && c<='9') || (c>='a' && c<='f'), "Blacklist sha80 must be lowercase hexadecimal");
				}
			}
			if(positionMode_.positional()){
				DDLFamilyProfiles.validatePositionalWidth(width_);
				DDLFamilyProfiles.validatePositionalLayout(partitions_,overlap_,geometry_,positionMode_);
			}
			geometry=geometry_; k=geometry.k(); width=width_; sketchSeed=sketchSeed_; positionMode=positionMode_; positional=positionMode_.positional();
			partitions=partitions_; overlap=overlap_; families=families_;
			blacklistSha80=blacklistSha80_;
			keyPolicy=keyPolicy_; keyEncodingVersion=keyEncodingVersion_; keyRawInputBinding=keyRawInputBinding_;
			consensusMemberBonus=consensusMemberBonus_;
			ddl.ProteinBlacklistPolicy.requireStored(geometry,blacklistSha80);
		}
		public Header withKeyPolicy(final DDLKeyTieBreak.PolicyContext context){
			require(context!=null, "Missing key-policy context");
			context.requireGeometry(width,geometry,sketchSeed,positionMode,partitions,overlap,blacklistSha80);
			return new Header(geometry,width,sketchSeed,positionMode,partitions,overlap,families,blacklistSha80,
				context.policy,context.encodingVersion,context.rawInputBinding,context.consensusMemberBonus);
		}
		public boolean sameKeyPolicy(final Header b){
			return b!=null && keyPolicy==b.keyPolicy && consensusMemberBonus==b.consensusMemberBonus &&
				java.util.Objects.equals(keyEncodingVersion,b.keyEncodingVersion) &&
				java.util.Objects.equals(keyRawInputBinding,b.keyRawInputBinding);
		}
		public boolean sameModelBinding(final Header b){return sameGeometry(b) && sameKeyPolicy(b);}
		public boolean sameGeometry(final Header b){
			return b!=null && geometry.equals(b.geometry) && width==b.width && sketchSeed==b.sketchSeed &&
				positionMode==b.positionMode && partitions==b.partitions && overlap==b.overlap &&
				java.util.Objects.equals(blacklistSha80,b.blacklistSha80);
		}
		public Header withFamilies(final int families_){
			return new Header(geometry,width,sketchSeed,positionMode,partitions,overlap,families_,blacklistSha80,
				keyPolicy,keyEncodingVersion,keyRawInputBinding,consensusMemberBonus);
		}
		public final int k, width, partitions, overlap, families;
		public final long sketchSeed;
		public final long consensusMemberBonus;
		public final boolean positional;
		public final ProteinPositionMode positionMode;
		public final ProteinKmerGeometry geometry;
		public final DDLKeyTieBreak.Policy keyPolicy;
		public final String blacklistSha80, keyEncodingVersion, keyRawInputBinding;
	}

	/** Raw model data; parsed arrays are owned, resident arrays are borrowed until import completes.
	 * Optional member descriptors are immutable resident metadata, never inferred from model payloads.
	 * Neither exported counts nor descriptor counts imply current assigned sequence IDs.
	 */
	public static final class Entry{
		private Entry(final int familyId_, final String name_, final int members_, final byte[] consensus_,
				final char[] sketch_, final float[] monomers_, final double[] scalars_, final int scalarMask_){
			this(familyId_,name_,name_,members_,consensus_,sketch_,monomers_,scalars_,scalarMask_,null);
		}
		private Entry(final int familyId_, final String name_, final String sketchIdentity_, final int members_, final byte[] consensus_,
				final char[] sketch_, final float[] monomers_, final double[] scalars_, final int scalarMask_,
				final MemberRankStatistics.Snapshot memberRanking_){
			validateMemberRanking(members_,memberRanking_);
			familyId=familyId_; name=name_; sketchIdentity=sketchIdentity_; members=members_; consensus=consensus_;
			sketch=sketch_; monomers=monomers_; scalars=scalars_; scalarMask=scalarMask_;
			memberRanking=memberRanking_;
			assert(consensus.length>0 && monomers.length==MONOMERS && scalars.length==SCALARS) :
				"Only validated complete records may construct a decoded entry";
		}
		public final int familyId, members, scalarMask;
		public final String name, sketchIdentity;
		public final byte[] consensus;
		public final char[] sketch;
		public final float[] monomers;
		public final double[] scalars;
		/** Null means unbound; it is never a measured zero or permission to substitute consensus statistics. */
		public final MemberRankStatistics.Snapshot memberRanking;
	}

	/** Header plus exactly its declared number of complete records. */
	public static final class Loaded{
		private Loaded(final Header header_, final Entry[] families_){
			header=header_; families=families_;
			assert(families.length==header.families) : "Declared and decoded family census must agree";
		}
		public final Header header;
		public final Entry[] families;
	}

	/** Rejects duplicate stable IDs or display names without boxing integer IDs. */
	private static final class Ids{
		void add(final int id, final String name){
			require(ids.get(id)<0 && names.add(name), "Duplicate cluster stable ID or name");
			ids.put(id,1);
		}
		private final IntHashMap ids=new IntHashMap(16);
		private final HashSet<String> names=new HashSet<String>();
	}

	private ClusterSketchCodec(){}
	public static final String SCHEMA="cluster_sketch_v2";
	public static final int LINES=5, MONOMERS=20, SCALARS=5;
	private static final int RANKING_FIELDS=4;
	private static final char[] HEX="0123456789ABCDEF".toCharArray();
	private static final byte[][] LABELS={
		"#consensus\tamino_acids".getBytes(StandardCharsets.US_ASCII),
		"#ddl\tA48_u16\tnonempty_mode".getBytes(StandardCharsets.US_ASCII),
		"#monomers\t20\tA48_f32".getBytes(StandardCharsets.US_ASCII),
		"#scalars\tmask_A48_u5\t5\tA48_f64\thydrophobicity\tnet_charge\taromaticity\tfraction_small\tpi_tail".getBytes(StandardCharsets.US_ASCII)
	};
}
