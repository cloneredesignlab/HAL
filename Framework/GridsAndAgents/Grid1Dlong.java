package Framework.GridsAndAgents;

import Framework.Tools.Internal.PDEequations;
import Framework.Util;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Created by Rafael on 10/24/2017.
 */
public class Grid1Dlong extends GridBase1D implements Serializable {
    long[] field;

    public Grid1Dlong(int xDim) {
        super(xDim, false);
        field = new long[this.xDim];
    }

    public Grid1Dlong(int xDim, boolean wrapX) {
        super(xDim, wrapX);
        field = new long[this.xDim];
    }

    /**
     * gets the current field value at the specified index
     */
    public long Get(int x) {
        return field[x];
    }

    /**
     * returns the complete field as an array
     */
    public long[] GetField() {
        return this.field;
    }

    /**
     * sets the current field value at the specified index
     */
    public void Set(int x, long val) {
        field[x] = val;
    }

    /**
     * multiplies the current field value at the specified index
     */
    public void Mul(int x, long val) {
        field[x] *= val;
    }

    /**
     * adds to the current field value at the specified index
     */
    public void Add(int x, long val) {
        field[x] += val;
    }

    /**
     * Bounds all values in the current field between min and max
     */
    public void BoundAll(long min, long max) {
        for (int i = 0; i < length; i++) {
            field[i] = Util.Bound(field[i], min, max);
        }
    }

    /**
     * sets all squares in current the field to the specified value
     */
    public void SetAll(long val) {
        Arrays.fill(field, val);
    }

    /**
     * adds specified value to all entries of the curr field
     */
    public void AddAll(long val) {
        for (int i = 0; i < length; i++) {
            field[i] += val;
        }
    }

    /**
     * adds specified value to all entries of the curr field
     */
    public void MulAll(long val) {
        for (int i = 0; i < length; i++) {
            field[i] *= val;
        }
    }

    /**
     * copies the array argument into the field
     */

    public void SetAll(long[] vals) {
        System.arraycopy(vals, 0, field, 0, length);
    }

    /**
     * returns the mean value of the grid
     */

    public long GetAvg() {
        long tot = 0;
        for (int i = 0; i < length; i++) {
            tot += field[i];
        }
        return tot / length;
    }

    /**
     * returns the max value in the grid
     */
    public long GetMax() {
        long max = Long.MIN_VALUE;
        for (int i = 0; i < length; i++) {
            max = Math.max(Get(i), max);
        }
        return max;
    }

    /**
     * returns the min value in the grid
     */
    public long GetMin() {
        long min = Long.MAX_VALUE;
        for (int i = 0; i < length; i++) {
            min = Math.min(Get(i), min);
        }
        return min;
    }
}