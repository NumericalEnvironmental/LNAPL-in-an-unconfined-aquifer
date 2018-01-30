# LNAPL-in-an-unconfined-aquifer

![Preview](https://numericalenvironmental.files.wordpress.com/2016/12/schematic.jpg?w=646&h=364)

A Julia-language script for modeling LNAPL floating on top of groundwater in a single-layer unconfined aquifer using an integral finite difference methodology. Specifically, the movement of immiscible fluids of contrasting density and effective hydraulic conductivity in porous media is simulated, including responses to extraction well operation and LNAPL bail-down testing. Input files include:

*nodes.txt – locations of volume element centroids, along with local initial conditions, hydraulic conductivity, and specific yield

*connects.txt – connection information for how the nodes are linked together in a model mesh (connection distances and interfaces)

*fluids.txt – fluid properties (just name and density for now)

*knobs.txt – basic model controls (time stepping constraints, etc.)

More info can be found here: https://numericalenvironmental.wordpress.com/2016/12/08/an-integrated-finite-difference-based-two-layer-lnapl-groundwater-model-written-in-julia/

I'd appreciate hearing back from you if you find the script useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
