package prot;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.zip.CRC32;

import fileIO.FileFormat;

/**
 * Stages a related output set and publishes it through same-filesystem atomic
 * moves. A filesystem-backed journal restores the complete old set after an
 * interrupted pre-commit publication, or retains the complete new set after a
 * commit marker is durable. Every output has its own stable lock, so overlapping
 * output sets cannot publish concurrently.
 *
 * <p>Staged files live in sibling directories and retain their exact final
 * basenames. This preserves BBTools compression selection and ZIP entry names.
 * The lock files intentionally remain on disk so all future processes lock the
 * same inode.
 */
final class OutputFileSetTransaction implements AutoCloseable {

	OutputFileSetTransaction(final String[] paths, final boolean overwrite_){
		this(paths,EMPTY_PATHS,overwrite_,MoveFault.NONE);
	}

	OutputFileSetTransaction(final String[] paths, final boolean overwrite_,
			final MoveFault fault_){
		this(paths,EMPTY_PATHS,overwrite_,fault_);
	}

	OutputFileSetTransaction(final String[] paths, final String[] protectedPaths,
			final boolean overwrite_){
		this(paths,protectedPaths,overwrite_,MoveFault.NONE);
	}

	OutputFileSetTransaction(final String[] paths, final String[] protectedPaths,
			final boolean overwrite_, final MoveFault fault_){
		overwrite=overwrite_; fault=fault_==null ? MoveFault.NONE : fault_;
		states=resolveStates(paths);
		validateStatePaths(states,protectedPaths);
		locks=LockSet.acquire(states);
		boolean initialized=false;
		try{
			validateExistingStateTypes(states);
			recoverLocked(states);
			assertNoOrphanState(states);
			hadOriginal=new boolean[states.length];
			for(int i=0; i<states.length; i++){
				final File target=states[i].target;
				if(target.exists() && (!overwrite || !target.isFile())){
					throw new IllegalArgumentException("Output exists: "+target);
				}
				hadOriginal[i]=target.exists();
			}
			writeJournal();
			for(final StatePath state : states){Files.createDirectory(state.stageDirectory.toPath());}
			forceParentDirectories(states);
			initialized=true;
		}catch(IOException e){throw new RuntimeException("Failed preparing output transaction",e);}
		finally{
			if(!initialized){
				try{recoverLocked(states);}catch(RuntimeException ignored){}
				locks.close();
			}
		}
	}

	/** Recovers an interrupted transaction before ordinary overwrite checks. */
	static void recover(final String[] paths, final String[] protectedPaths){
		final StatePath[] states=resolveStates(paths);
		validateStatePaths(states,protectedPaths);
		final LockSet locks=LockSet.acquire(states);
		try{validateExistingStateTypes(states); recoverLocked(states);}finally{locks.close();}
	}

	/** Test and maintenance shorthand when no input paths need protection. */
	static void recover(final String[] paths){recover(paths,EMPTY_PATHS);}

	String temporaryPath(final int index){return states[index].temporary.getPath();}

	/** Publishes the entire staged set or restores the entire preceding set. */
	void commit(){
		validateExistingStateTypes(states);
		for(final StatePath state : states){
			if(!state.temporary.isFile()){
				throw new RuntimeException("Missing staged output: "+state.temporary);
			}
		}
		try{
			for(final StatePath state : states){
				forceFile(state.temporary); forceDirectory(state.stageDirectory);
			}
			for(int i=0; i<states.length; i++){
				final StatePath state=states[i];
				if(hadOriginal[i]){
					if(!state.target.isFile()){
						throw new IOException("Original output disappeared during write: "+state.target);
					}
					fault.beforeMove("backup",i);
					moveAtomic(state.target,state.backup);
					fault.afterMove("backup",i);
				}else if(state.target.exists()){
					throw new IOException("Output appeared during write: "+state.target);
				}
			}
			forceParentDirectories(states);
			for(int i=0; i<states.length; i++){
				fault.beforeMove("publish",i);
				moveAtomic(states[i].temporary,states[i].target);
				fault.afterMove("publish",i);
			}
			for(final StatePath state : states){forceDirectory(state.stageDirectory);}
			forceParentDirectories(states);
			writeCommitMarker(states[0]);
			committed=true;
		}catch(Throwable failure){
			if(commitMarkerVisible(states[0])){
				committed=true;
				System.err.println("WARNING: commit marker exists after publication exception; "+
					"retaining the complete new output set: "+failure.getMessage());
				finishCommittedCleanup();
				return;
			}
			try{recoverLocked(states);}catch(Throwable recoveryFailure){failure.addSuppressed(recoveryFailure);}
			if(failure instanceof Error){throw (Error)failure;}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw new RuntimeException("Output transaction failed",failure);
		}
		finishCommittedCleanup();
	}

	@Override
	public void close(){
		RuntimeException failure=null;
		if(!committed){
			try{recoverLocked(states);}catch(RuntimeException e){failure=e;}
		}
		try{locks.close();}
		catch(RuntimeException e){
			if(failure==null){failure=e;}else{failure.addSuppressed(e);}
		}
		if(failure!=null){throw failure;}
	}

	private void finishCommittedCleanup(){
		try{
			fault.beforeCleanup();
			recoverLocked(states);
		}catch(Exception cleanupFailure){
			System.err.println("WARNING: outputs are committed; deferred transaction cleanup: "+
				cleanupFailure.getMessage());
		}
	}

	private void writeJournal() throws IOException{
		final StatePath primary=states[0];
		Files.deleteIfExists(primary.journalPending.toPath());
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final DataOutputStream payload=new DataOutputStream(buffer);
		payload.writeInt(JOURNAL_MAGIC); payload.writeInt(JOURNAL_VERSION);
		payload.writeInt(states.length);
		for(int i=0; i<states.length; i++){
			payload.writeBoolean(hadOriginal[i]);
			payload.writeUTF(states[i].target.getPath());
		}
		payload.flush();
		final byte[] bytes=buffer.toByteArray();
		final CRC32 crc=new CRC32(); crc.update(bytes);
		final FileOutputStream fos=new FileOutputStream(primary.journalPending);
		try{
			final DataOutputStream out=new DataOutputStream(new BufferedOutputStream(fos));
			out.write(bytes); out.writeLong(crc.getValue());
			out.flush(); fos.getFD().sync(); out.close();
		}finally{try{fos.close();}catch(IOException ignored){}}
		moveAtomic(primary.journalPending,primary.journal);
		forceDirectory(primary.journal.getParentFile());
	}

	private static void writeCommitMarker(final StatePath primary) throws IOException{
		Files.deleteIfExists(primary.markerPending.toPath());
		final FileOutputStream fos=new FileOutputStream(primary.markerPending);
		try{
			fos.write(COMMITTED_BYTES); fos.flush(); fos.getFD().sync();
		}finally{fos.close();}
		moveAtomic(primary.markerPending,primary.committedMarker);
		forceDirectory(primary.committedMarker.getParentFile());
	}

	private static void recoverLocked(final StatePath[] requested){
		final StatePath primary=requested[0];
		try{
			validateExistingStateTypes(requested);
			Files.deleteIfExists(primary.journalPending.toPath());
			Files.deleteIfExists(primary.markerPending.toPath());
			if(!primary.journal.exists()){
				Files.deleteIfExists(primary.committedMarker.toPath());
				return;
			}
			final JournalRecord record=readJournal(primary.journal,requested);
			final boolean committed=validCommitMarker(primary.committedMarker);
			if(committed){
				for(final StatePath state : requested){
					if(!state.target.isFile()){
						throw new IOException("Committed transaction is missing final output: "+state.target);
					}
				}
				for(final StatePath state : requested){
					Files.deleteIfExists(state.temporary.toPath());
					deleteEmptyDirectory(state.stageDirectory);
					Files.deleteIfExists(state.backup.toPath());
				}
				forceParentDirectories(requested);
				Files.delete(primary.journal.toPath());
				forceDirectory(primary.journal.getParentFile());
				Files.deleteIfExists(primary.committedMarker.toPath());
				forceDirectory(primary.committedMarker.getParentFile());
				return;
			}
			Files.deleteIfExists(primary.committedMarker.toPath());
			for(int i=requested.length-1; i>=0; i--){
				final StatePath state=requested[i];
				if(record.hadOriginal[i]){
					if(state.backup.exists()){
						Files.deleteIfExists(state.target.toPath());
						moveAtomic(state.backup,state.target);
					}
				}else if(state.stageDirectory.exists()){
					if(state.backup.exists()){
						throw new IOException("Unexpected backup for originally absent output: "+state.backup);
					}
					Files.deleteIfExists(state.target.toPath());
				}
				Files.deleteIfExists(state.temporary.toPath());
				deleteEmptyDirectory(state.stageDirectory);
			}
			forceParentDirectories(requested);
			Files.delete(primary.journal.toPath());
			forceDirectory(primary.journal.getParentFile());
		}catch(IOException e){
			throw new RuntimeException("Failed recovering output transaction "+primary.journal,e);
		}
	}

	private static JournalRecord readJournal(final File journal,
			final StatePath[] requested) throws IOException{
		final byte[] bytes=Files.readAllBytes(journal.toPath());
		if(bytes.length<20){throw new IOException("Truncated output transaction journal: "+journal);}
		final int payloadLength=bytes.length-8;
		final DataInputStream checksumIn=new DataInputStream(
			new ByteArrayInputStream(bytes,payloadLength,8));
		final long expectedChecksum=checksumIn.readLong(); checksumIn.close();
		final CRC32 crc=new CRC32(); crc.update(bytes,0,payloadLength);
		if(crc.getValue()!=expectedChecksum){
			throw new IOException("Output transaction journal checksum mismatch: "+journal);
		}
		final DataInputStream in=new DataInputStream(
			new BufferedInputStream(new ByteArrayInputStream(bytes,0,payloadLength)));
		try{
			if(in.readInt()!=JOURNAL_MAGIC || in.readInt()!=JOURNAL_VERSION){
				throw new IOException("Invalid output transaction journal header: "+journal);
			}
			final int count=in.readInt();
			if(count!=requested.length){
				throw new IOException("Journal output count does not match requested set: "+count);
			}
			final boolean[] existed=new boolean[count];
			for(int i=0; i<count; i++){
				existed[i]=in.readBoolean();
				final String path=in.readUTF();
				if(!path.equals(requested[i].target.getPath())){
					throw new IOException("Journal output does not match requested set at index "+i);
				}
			}
			if(in.read()!=-1){throw new IOException("Trailing bytes in output transaction journal: "+journal);}
			return new JournalRecord(existed);
		}finally{in.close();}
	}

	private static StatePath[] resolveStates(final String[] paths){
		if(paths==null || paths.length<1){throw new IllegalArgumentException("No output paths.");}
		final StatePath[] states=new StatePath[paths.length];
		final HashSet<String> seen=new HashSet<String>();
		try{
			for(int i=0; i<paths.length; i++){
				final String path=paths[i];
				if(path==null || path.length()==0 || FileFormat.isStdio(path)){
					throw new IllegalArgumentException("Transactional output must be a file: "+path);
				}
				final File target=new File(path).getAbsoluteFile();
				final File parent=target.getParentFile();
				if(parent==null || !parent.isDirectory()){
					throw new IllegalArgumentException("Output directory does not exist: "+path);
				}
				if(!seen.add(target.getCanonicalPath())){
					throw new IllegalArgumentException("Duplicate output: "+path);
				}
				states[i]=new StatePath(target);
			}
		}catch(IOException e){throw new RuntimeException("Failed resolving output paths",e);}
		return states;
	}

	private static void validateStatePaths(final StatePath[] states,
			final String[] protectedPaths){
		final ArrayList<File> stateFiles=new ArrayList<File>();
		for(final StatePath state : states){state.addFiles(stateFiles);}
		final HashSet<String> stateCanonical=new HashSet<String>();
		try{
			for(final File file : stateFiles){
				final String canonical=file.getCanonicalPath();
				if(!stateCanonical.add(canonical)){
					throw new IllegalArgumentException("Transaction state path collision: "+file);
				}
			}
			for(int i=0; i<stateFiles.size(); i++){
				final File a=stateFiles.get(i);
				if(!a.exists()){continue;}
				for(int j=i+1; j<stateFiles.size(); j++){
					final File b=stateFiles.get(j);
					if(b.exists() && Files.isSameFile(a.toPath(),b.toPath())){
						throw new IllegalArgumentException(
							"Transaction state paths are hard-linked: "+a+" and "+b);
					}
				}
			}
			if(protectedPaths!=null){
				for(final String path : protectedPaths){
					if(path==null){continue;}
					final File protectedFile=new File(path);
					final String canonical=protectedFile.getCanonicalPath();
					if(stateCanonical.contains(canonical)){
						throw new IllegalArgumentException(
							"Input collides with transaction state path: "+path);
					}
					if(protectedFile.exists()){
						for(final File stateFile : stateFiles){
							if(stateFile.exists() && Files.isSameFile(protectedFile.toPath(),stateFile.toPath())){
								throw new IllegalArgumentException(
									"Input hard-links transaction state path: "+path);
							}
						}
					}
				}
			}
		}catch(IOException e){throw new RuntimeException("Failed resolving transaction state paths",e);}
	}

	private static void assertNoOrphanState(final StatePath[] states) throws IOException{
		for(int i=0; i<states.length; i++){
			final StatePath state=states[i];
			if(state.stageDirectory.exists() || state.backup.exists() ||
					(i>0 && (state.journal.exists() || state.committedMarker.exists() ||
					state.journalPending.exists() || state.markerPending.exists()))){
				throw new IOException("Orphan state belongs to another output transaction: "+state.target);
			}
		}
	}

	private static void validateExistingStateTypes(final StatePath[] states){
		for(final StatePath state : states){
			validateDirectory(state.stageDirectory);
			validateRegularFile(state.temporary);
			validateRegularFile(state.backup);
			validateRegularFile(state.lockFile);
			validateRegularFile(state.journal);
			validateRegularFile(state.journalPending);
			validateRegularFile(state.committedMarker);
			validateRegularFile(state.markerPending);
		}
	}

	private static void validateDirectory(final File file){
		if(file.exists() && (Files.isSymbolicLink(file.toPath()) || !file.isDirectory())){
			throw new IllegalArgumentException("Unsafe transaction directory: "+file);
		}
	}

	private static void validateRegularFile(final File file){
		if(file.exists() && (Files.isSymbolicLink(file.toPath()) || !file.isFile())){
			throw new IllegalArgumentException("Unsafe transaction file: "+file);
		}
	}

	private static void forceFile(final File file) throws IOException{
		final FileChannel channel=FileChannel.open(file.toPath(),StandardOpenOption.WRITE);
		try{channel.force(true);}finally{channel.close();}
	}

	private static void forceParentDirectories(final StatePath[] states) throws IOException{
		final HashSet<String> seen=new HashSet<String>();
		for(final StatePath state : states){
			final File parent=state.target.getParentFile();
			if(seen.add(parent.getCanonicalPath())){forceDirectory(parent);}
		}
	}

	private static void forceDirectory(final File directory) throws IOException{
		final FileChannel channel=FileChannel.open(directory.toPath(),StandardOpenOption.READ);
		try{channel.force(true);}finally{channel.close();}
	}

	private static void moveAtomic(final File from, final File to) throws IOException{
		Files.move(from.toPath(),to.toPath(),StandardCopyOption.ATOMIC_MOVE);
	}

	private static void deleteEmptyDirectory(final File directory) throws IOException{
		if(directory.exists()){Files.delete(directory.toPath());}
	}

	private static boolean validCommitMarker(final File marker) throws IOException{
		if(!marker.isFile()){return false;}
		return Arrays.equals(Files.readAllBytes(marker.toPath()),COMMITTED_BYTES);
	}

	private static boolean commitMarkerVisible(final StatePath primary){
		try{return validCommitMarker(primary.committedMarker);}catch(IOException e){return false;}
	}

	/** Test-only fault boundary; production uses {@link #NONE}. */
	interface MoveFault{
		void beforeMove(String phase, int index) throws IOException;
		void afterMove(String phase, int index) throws IOException;
		void beforeCleanup() throws IOException;
		MoveFault NONE=new MoveFault(){
			@Override public void beforeMove(final String phase, final int index){}
			@Override public void afterMove(final String phase, final int index){}
			@Override public void beforeCleanup(){}
		};
	}

	private static final class StatePath{
		StatePath(final File target_){
			target=target_;
			final File parent=target.getParentFile();
			final String hidden="."+target.getName()+".bbtxn";
			stageDirectory=new File(parent,hidden+".stage");
			temporary=new File(stageDirectory,target.getName());
			backup=new File(parent,hidden+".backup");
			lockFile=new File(parent,hidden+".lock");
			journal=new File(parent,hidden+".journal");
			journalPending=new File(parent,hidden+".journal.pending");
			committedMarker=new File(parent,hidden+".committed");
			markerPending=new File(parent,hidden+".committed.pending");
		}
		void addFiles(final ArrayList<File> files){
			files.add(target); files.add(stageDirectory); files.add(temporary); files.add(backup);
			files.add(lockFile); files.add(journal); files.add(journalPending);
			files.add(committedMarker); files.add(markerPending);
		}
		final File target,stageDirectory,temporary,backup,lockFile;
		final File journal,journalPending,committedMarker,markerPending;
	}

	private static final class JournalRecord{
		JournalRecord(final boolean[] existed){hadOriginal=existed;}
		final boolean[] hadOriginal;
	}

	private static final class LockSet implements AutoCloseable{
		static LockSet acquire(final StatePath[] states){
			final File[] files=new File[states.length];
			for(int i=0; i<states.length; i++){files[i]=states[i].lockFile;}
			Arrays.sort(files,new Comparator<File>(){
				@Override
				public int compare(final File a, final File b){return a.getPath().compareTo(b.getPath());}
			});
			final LockSet set=new LockSet(files);
			try{set.acquireAll(); return set;}
			catch(RuntimeException e){set.close(); throw e;}
		}

		LockSet(final File[] files){
			lockFiles=files; channels=new FileChannel[files.length]; locks=new FileLock[files.length];
		}

		private void acquireAll(){
			try{
				for(int i=0; i<lockFiles.length; i++){
					channels[i]=FileChannel.open(lockFiles[i].toPath(),StandardOpenOption.CREATE,
						StandardOpenOption.WRITE);
					try{locks[i]=channels[i].tryLock();}
					catch(OverlappingFileLockException e){
						throw new IOException("Output transaction is active: "+lockFiles[i],e);
					}
					if(locks[i]==null){throw new IOException("Output transaction is active: "+lockFiles[i]);}
				}
			}catch(IOException e){throw new RuntimeException("Failed locking output transaction",e);}
		}

		@Override
		public void close(){
			RuntimeException failure=null;
			for(int i=locks.length-1; i>=0; i--){
				try{if(locks[i]!=null){locks[i].release();}}
				catch(IOException e){
					if(failure==null){failure=new RuntimeException("Failed releasing transaction lock",e);}
					else{failure.addSuppressed(e);}
				}
				try{if(channels[i]!=null){channels[i].close();}}
				catch(IOException e){
					if(failure==null){failure=new RuntimeException("Failed closing transaction lock",e);}
					else{failure.addSuppressed(e);}
				}
			}
			if(failure!=null){throw failure;}
		}

		private final File[] lockFiles;
		private final FileChannel[] channels;
		private final FileLock[] locks;
	}

	private final StatePath[] states;
	private boolean[] hadOriginal;
	private final boolean overwrite;
	private final MoveFault fault;
	private final LockSet locks;
	private boolean committed;
	private static final int JOURNAL_MAGIC=0x42425458;
	private static final int JOURNAL_VERSION=3;
	private static final byte[] COMMITTED_BYTES="BBTOOLS_OUTPUT_COMMITTED_V2\n".
		getBytes(StandardCharsets.US_ASCII);
	private static final String[] EMPTY_PATHS=new String[0];
}
