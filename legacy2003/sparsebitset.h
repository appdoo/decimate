/******************************************************************************
 * 
 * SparseBitSet
 * 
 *
 * A sparse Bit-Set implementation
 * specially for memory saving storage of indes lists
 *
 *****************************************************************************/
#define BLOCKBASE 6
#define BLOCKSIZE (1 << BLOCKBASE)

class SparseBitSet
{
    int  *offs;
    long *bits;
    int size;
public:
    /**************************************************************************
     * Creates an empty bitset.
     *************************************************************************/
    public SparseBitSet() 
    {
   	bits = new long[4];
	   offs = new int [4];
	   size = 0;
    }

    private void new_block(int bnum) 
    {
   	new_block(bsearch(bnum), bnum);
    }

   private void new_block(int idx, int bnum) 
   {
	   if (size==bits.length) { // resize
	      long[] nbits = new long[size*3];
	      int [] noffs = new int [size*3];
	      System.arraycopy(bits, 0, nbits, 0, size);
	      System.arraycopy(offs, 0, noffs, 0, size);
	      bits = nbits;
	      offs = noffs;
	   }
   	insert_block(idx, bnum);
   }

   private void insert_block(int idx, int bnum) 
   {
	   CUtility.ASSERT(idx<=size);
	   CUtility.ASSERT(idx==size || offs[idx]!=bnum);
	   System.arraycopy(bits, idx, bits, idx+1, size-idx);
	   System.arraycopy(offs, idx, offs, idx+1, size-idx);
	   offs[idx]=bnum;
	   bits[idx]=0; //clear them bits.
	   size++;
   }

   private int bsearch(int bnum) 
   {
	   int l=0, r=size; // search interval is [l, r)
      while (l<r) {
	      int p = (l+r)/2;
	      if (bnum<offs[p]) r=p;
	      else if (bnum>offs[p]) l=p+1;
	      else return p;
   	}
	   CUtility.ASSERT(l==r);
	   return l; // index at which the bnum *should* be, if it's not.
   }

   /**
    * Sets a bit.
    * @param bit the bit to be set
    */
   public void set(int bit) 
   {
	   int bnum = bit >> LG_BITS;
	   int idx  = bsearch(bnum);
	   if (idx >= size || offs[idx]!=bnum)
	      new_block(idx, bnum);
	   bits[idx] |= (1L << (bit & BITS_M1) );
   }

   /**
    * Clears a bit.
    * @param bit the bit to be cleared
    */
   public void clear(int bit) 
   {
	   int bnum = bit >> LG_BITS;
	   int idx  = bsearch(bnum);
	   if (idx >= size || offs[idx]!=bnum)
	      new_block(idx, bnum);
	   bits[idx] &= ~(1L << (bit & BITS_M1) );
   }

   /**
    * Clears all bits.
    */
   public void clearAll() 
   {
	   size = 0;
   }

   /**
    * Gets a bit.
    * @param bit the bit to be gotten
    */
   public boolean get(int bit) 
   {
	   int bnum = bit >> LG_BITS;
	   int idx  = bsearch(bnum);
	   if (idx >= size || offs[idx]!=bnum)
	      return false;
	   return 0 != ( bits[idx] & (1L << (bit & BITS_M1) ) );
   }

   /**
    * Logically ANDs this bit set with the specified set of bits.
    * @param set the bit set to be ANDed with
    */
   public void and(SparseBitSet set) 
   {
	   binop(this, set, AND);
   }

   /**
    * Logically ORs this bit set with the specified set of bits.
    * @param set the bit set to be ORed with
    */
   void or(SparseBitSet set) 
   {
	   binop(this, set, OR);
   }

   /**
    * Logically XORs this bit set with the specified set of bits.
    * @param set the bit set to be XORed with
    */
   void xor(SparseBitSet set) 
   {
	   binop(this, set, XOR);
   }

    // BINARY OPERATION MACHINERY
    private static interface BinOp {
	public long op(long a, long b);
    }
    private static final BinOp AND = new BinOp() {
	public final long op(long a, long b) { return a & b; }
    };
    private static final BinOp OR = new BinOp() {
	public final long op(long a, long b) { return a | b; }
    };
    private static final BinOp XOR = new BinOp() {
	public final long op(long a, long b) { return a ^ b; }
    };
    private static final void binop(SparseBitSet a, SparseBitSet b, BinOp op) {
	int  nsize = a.size + b.size;
	long[] nbits; 
	int [] noffs;
	int a_zero, a_size;
	// be very clever and avoid allocating more memory if we can.
	if (a.bits.length < nsize) { // oh well, have to make working space.
	    nbits = new long[nsize];
	    noffs = new int [nsize];
	    a_zero  = 0; a_size = a.size;
	} else { // reduce, reuse, recycle!
	    nbits = a.bits;
	    noffs = a.offs;
	    a_zero = a.bits.length - a.size; a_size = a.bits.length;
	    System.arraycopy(a.bits, 0, a.bits, a_zero, a.size);
	    System.arraycopy(a.offs, 0, a.offs, a_zero, a.size);
	}
	// ok, crunch through and binop those sets!
	nsize = 0;
	for (int i=a_zero, j=0; i<a_size || j<b.size; ) {
	    long nb; int no;
	    if (i<a_size && (j>=b.size || a.offs[i] < b.offs[j])) {
		nb = op.op(a.bits[i], 0);
		no = a.offs[i];
		i++;
	    } else if (j<b.size && (i>=a_size || a.offs[i] > b.offs[j])) {
		nb = op.op(0, b.bits[j]);
		no = b.offs[j];
		j++;
	    } else { // equal keys; merge.
		nb = op.op(a.bits[i], b.bits[j]);
		no = a.offs[i];
		i++; j++;
	    }
	    if (nb!=0) {
		nbits[nsize] = nb;
		noffs[nsize] = no;
		nsize++;
	    }
	}
	a.bits = nbits;
	a.offs = noffs;
	a.size = nsize;
    }

   /**
    * Calculates and returns the set's size
    */
   int size() 
   {
      return (size==0)?0:((1+offs[size-1]) << LG_BITS);
   }
}

